"""
Test suite for beamy.analysis module.
Tests FEM solver, reaction calculations, and analysis results.
"""

import numpy as np
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from beamy.analysis import (
    solve_x_reactions,
    solve_transverse_reactions,
    get_all_loads,
    Result,
    AnalysisResult,
    LoadedBeam,
)
from beamy.analysis.analysis import (
    _solve_fem_1d,
    _accumulate_loads,
    _hermite_displacement,
    _linear_interpolation,
)
from beamy.setup import Beam1D, Material, Support, LoadCase, PointForce, Moment
from beamy.section import Section


def test_result_class():
    """Test Result wrapper class."""
    print("Testing Result class...")
    x = np.array([0.0, 1.0, 2.0, 3.0])
    values = np.array([10.0, 20.0, 15.0, 5.0])
    r = Result(x, values)
    
    # Test iteration
    assert list(r) == [(0.0, 10.0), (1.0, 20.0), (2.0, 15.0), (3.0, 5.0)]
    
    # Test properties
    assert r.max == 20.0
    assert r.min == 5.0
    assert abs(r.mean - 12.5) < 1e-9
    assert r.range == 15.0
    
    # Test indexing
    assert r[0] == (0.0, 10.0)
    assert r[1] == (1.0, 20.0)
    
    # Test interpolation
    assert abs(r.at(0.5) - 15.0) < 1e-9  # Linear interpolation between 10 and 20
    assert abs(r.at(1.5) - 17.5) < 1e-9  # Linear interpolation between 20 and 15
    
    print("  [PASS] Result class tests passed")


def test_analysis_result():
    """Test AnalysisResult wrapper."""
    print("Testing AnalysisResult class...")
    x = np.array([0.0, 1.0, 2.0])
    action = Result(x, np.array([10.0, 20.0, 15.0]))
    stress = Result(x, np.array([100.0, 200.0, 150.0]))
    displacement = Result(x, np.array([0.0, 0.1, 0.2]))
    
    ar = AnalysisResult(action, stress, displacement)
    
    assert ar.action is action
    assert ar.stress is stress
    assert ar.displacement is displacement
    
    print("  [PASS] AnalysisResult class tests passed")


def test_accumulate_loads():
    """Test load accumulation function."""
    print("Testing _accumulate_loads...")
    points = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    loads = [(1.0, 100.0), (3.0, -50.0)]  # Force at x=1, force at x=3
    moment_loads = [(2.0, 200.0)]  # Moment at x=2
    
    force_vals, moment_vals = _accumulate_loads(points, loads, moment_loads)
    
    # At x=0: no loads
    assert force_vals[0] == 0.0
    assert moment_vals[0] == 0.0
    
    # At x=1: one force (100)
    assert force_vals[1] == 100.0
    assert abs(moment_vals[1] - 0.0) < 1e-9  # No moment yet
    
    # At x=2: one force (100), one moment (200) + moment from force
    assert force_vals[2] == 100.0
    assert abs(moment_vals[2] - (200.0 + 100.0 * 1.0)) < 1e-9  # Moment + force * arm
    
    # At x=3: two forces (100, -50), moment (200) + moments from forces
    assert force_vals[3] == 50.0  # 100 - 50
    expected_moment = 200.0 + 100.0 * 2.0 + (-50.0) * 0.0
    assert abs(moment_vals[3] - expected_moment) < 1e-9
    
    # Test without moment loads (forces still contribute to moments via moment arm)
    force_vals2, moment_vals2 = _accumulate_loads(points, loads, moment_loads=[])
    assert force_vals2[2] == 100.0
    # At x=2, force at x=1 contributes moment: 100 * (2-1) = 100
    assert abs(moment_vals2[2] - 100.0 * 1.0) < 1e-9
    
    print("  [PASS] _accumulate_loads tests passed")


def test_linear_interpolation():
    """Test linear interpolation function."""
    print("Testing _linear_interpolation...")
    points = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
    x_nodes = np.array([0.0, 1.0, 2.0])
    d_vec = np.array([0.0, 10.0, 20.0])  # Displacements at nodes
    
    result = _linear_interpolation(points, x_nodes, d_vec)
    
    assert abs(result[0] - 0.0) < 1e-9
    assert abs(result[1] - 5.0) < 1e-9  # Midpoint between 0 and 10
    assert abs(result[2] - 10.0) < 1e-9
    assert abs(result[3] - 15.0) < 1e-9  # Midpoint between 10 and 20
    assert abs(result[4] - 20.0) < 1e-9
    
    print("  [PASS] _linear_interpolation tests passed")


def test_hermite_displacement():
    """Test Hermite interpolation function."""
    print("Testing _hermite_displacement...")
    points = np.array([0.0, 0.5, 1.0])
    x_nodes = np.array([0.0, 1.0])
    # d_vec: [w0, theta0, w1, theta1]
    d_vec = np.array([0.0, 0.0, 1.0, 0.0])  # Simple case: w(0)=0, w(1)=1, no rotation
    
    result = _hermite_displacement(points, x_nodes, d_vec)
    
    assert abs(result[0] - 0.0) < 1e-6
    assert abs(result[2] - 1.0) < 1e-6
    # At midpoint, should be approximately 0.5 for linear displacement
    assert abs(result[1] - 0.5) < 1e-3
    
    print("  [PASS] _hermite_displacement tests passed")


def test_solve_fem_1d_simple():
    """Test the generic FEM solver with a simple case."""
    print("Testing _solve_fem_1d (simple case)...")
    
    # Simple 2-node axial problem: fixed at x=0, force at x=1
    nodes = [
        Support(x=0.0, type="100000"),  # Fixed in x
        Support(x=1.0, type="000000"),  # Free
    ]
    
    def k_linear(L):
        k = 1.0 / L
        return k * np.array([[1.0, -1.0], [-1.0, 1.0]])
    
    def f_global(nodes, ndof):
        f = np.zeros(ndof)
        f[1] = 100.0  # Force at node 1
        return f
    
    def is_fixed(support, k):
        return support.type[0] == "1"  # Fixed in x-direction
    
    d, r = _solve_fem_1d(nodes, dof_per_node=1, k_local_fn=k_linear, 
                         f_global_fn=f_global, is_fixed_fn=is_fixed)
    
    # Node 0 should be fixed (d[0] = 0)
    assert abs(d[0]) < 1e-9
    # Node 1 should have displacement
    assert d[1] > 0
    
    # Reaction at node 0 should balance the force
    assert abs(r[0] + 100.0) < 1e-6  # Reaction should be -100
    
    print("  [PASS] _solve_fem_1d simple tests passed")


def test_solve_x_reactions():
    """Test axial and torsional reaction solver."""
    print("Testing solve_x_reactions...")
    
    # Create a simple beam with two supports
    # Need at least one support in x, y, z, and rotation about x
    supports = [
        Support(x=0.0, type="111100"),  # Fixed in x, y, z, rotation about x
        Support(x=5.0, type="011100"),  # Free in x, fixed others
    ]
    
    # Create load case with axial force
    loads = LoadCase(name="Test")
    loads.add_point_force(PointForce(
        point=np.array([2.5, 0.0, 0.0]),
        force=np.array([100.0, 0.0, 0.0])  # 100N in x-direction
    ))
    
    d_x, d_rx = solve_x_reactions(supports, loads)
    
    # Check that reactions were computed
    assert abs(supports[0].reactions["Fx"] + 100.0) < 1e-6  # Reaction should balance
    assert abs(supports[1].reactions["Fx"]) < 1e-6  # Free end, no reaction
    
    # Check displacements
    assert len(d_x) == 2
    assert len(d_rx) == 2
    
    print("  [PASS] solve_x_reactions tests passed")


def test_solve_transverse_reactions():
    """Test transverse reaction solver."""
    print("Testing solve_transverse_reactions...")
    
    # Create beam with material and section
    material = Material(name="Steel", E=200e9, G=80e9)
    section = Section(
        name="Rect",
        A=0.01,  # 10cm x 10cm
        Iy=1e-6,  # Iy for z-bending
        Iz=1e-6,  # Iz for y-bending
        J=2e-6,
        y_max=0.05,
        z_max=0.05,
    )
    
    supports = [
        Support(x=0.0, type="111100"),  # Fixed translations and rotation about x
        Support(x=5.0, type="111100"),  # Fixed translations and rotation about x
    ]
    
    beam = Beam1D(L=5.0, material=material, section=section, supports=supports)
    
    # Create load case with transverse force
    loads = LoadCase(name="Test")
    loads.add_point_force(PointForce(
        point=np.array([2.5, 0.0, 0.0]),
        force=np.array([0.0, 0.0, -1000.0])  # 1000N downward (z-direction)
    ))
    
    d = solve_transverse_reactions(beam, loads, axis="z")
    
    # Check that reactions were computed
    # Total reaction should balance the applied force
    total_reaction = supports[0].reactions["Fz"] + supports[1].reactions["Fz"]
    assert abs(total_reaction - 1000.0) < 1e-6
    
    # Check displacement vector length (2 DOFs per node)
    assert len(d) == 4  # 2 nodes * 2 DOFs
    
    print("  [PASS] solve_transverse_reactions tests passed")


def test_get_all_loads():
    """Test load aggregation function."""
    print("Testing get_all_loads...")
    
    material = Material(name="Steel", E=200e9, G=80e9)
    section = Section(
        name="Rect",
        A=0.01,
        Iy=1e-6,
        Iz=1e-6,
        J=2e-6,
        y_max=0.05,
        z_max=0.05,
    )
    
    supports = [
        Support(x=0.0, type="111111"),  # Fully fixed
        Support(x=5.0, type="111111"),  # Fully fixed
    ]
    
    # Set some reactions manually
    supports[0].reactions["Fx"] = 50.0
    supports[1].reactions["Fz"] = 100.0
    
    beam = Beam1D(L=5.0, material=material, section=section, supports=supports)
    
    loads = LoadCase(name="Test")
    loads.add_point_force(PointForce(
        point=np.array([2.5, 0.0, 0.0]),
        force=np.array([100.0, 0.0, -200.0])
    ))
    
    all_loads = get_all_loads(loads, beam)
    
    # Check that applied loads are included
    fx_loads = [v for x, t, v in all_loads if t == "Fx"]
    assert 100.0 in fx_loads
    
    fz_loads = [v for x, t, v in all_loads if t == "Fz"]
    assert -200.0 in fz_loads
    
    # Check that reactions are included
    rx_loads = [v for x, t, v in all_loads if t == "Rx"]
    assert 50.0 in rx_loads
    
    rz_loads = [v for x, t, v in all_loads if t == "Rz"]
    assert 100.0 in rz_loads
    
    # Check sorting
    xs = [x for x, t, v in all_loads]
    assert xs == sorted(xs)
    
    print("  [PASS] get_all_loads tests passed")


def test_loaded_beam_basic():
    """Test LoadedBeam basic functionality."""
    print("Testing LoadedBeam (basic)...")
    
    material = Material(name="Steel", E=200e9, G=80e9)
    section = Section(
        name="Rect",
        A=0.01,
        Iy=1e-6,
        Iz=1e-6,
        J=2e-6,
        y_max=0.05,
        z_max=0.05,
    )
    
    supports = [
        Support(x=0.0, type="111111"),  # Fully fixed
        Support(x=5.0, type="111100"),  # Fixed translations and rotation about x
    ]
    
    beam = Beam1D(L=5.0, material=material, section=section, supports=supports)
    
    loads = LoadCase(name="Test")
    loads.add_point_force(PointForce(
        point=np.array([2.5, 0.0, 0.0]),
        force=np.array([100.0, 0.0, -1000.0])
    ))
    
    loaded_beam = LoadedBeam(beam, loads)
    
    # Check that all_loads was populated
    assert len(loaded_beam.all_loads) > 0
    
    # Test analysis methods
    shear_result = loaded_beam.shear("z", points=10)
    assert isinstance(shear_result, AnalysisResult)
    assert len(list(shear_result.action)) == 10
    
    bending_result = loaded_beam.bending("z", points=10)
    assert isinstance(bending_result, AnalysisResult)
    
    axial_result = loaded_beam.axial(points=10)
    assert isinstance(axial_result, AnalysisResult)
    
    torsion_result = loaded_beam.torsion(points=10)
    assert isinstance(torsion_result, AnalysisResult)
    
    deflection = loaded_beam.deflection("z", points=10)
    assert isinstance(deflection, Result)
    
    print("  [PASS] LoadedBeam basic tests passed")


def test_loaded_beam_equilibrium():
    """Test that LoadedBeam maintains equilibrium."""
    print("Testing LoadedBeam equilibrium...")
    
    material = Material(name="Steel", E=200e9, G=80e9)
    section = Section(
        name="Rect",
        A=0.01,
        Iy=1e-6,
        Iz=1e-6,
        J=2e-6,
        y_max=0.05,
        z_max=0.05,
    )
    
    supports = [
        Support(x=0.0, type="111111"),  # Fully fixed
        Support(x=5.0, type="111100"),  # Fixed translations and rotation about x
    ]
    
    beam = Beam1D(L=5.0, material=material, section=section, supports=supports)
    
    loads = LoadCase(name="Test")
    # Apply 1000N downward at midspan
    loads.add_point_force(PointForce(
        point=np.array([2.5, 0.0, 0.0]),
        force=np.array([0.0, 0.0, -1000.0])
    ))
    
    loaded_beam = LoadedBeam(beam, loads)
    
    # Check force equilibrium in z-direction
    total_applied = sum(v for x, t, v in loaded_beam.all_loads if t == "Fz")
    total_reactions = sum(v for x, t, v in loaded_beam.all_loads if t == "Rz")
    
    # Applied + reactions should sum to zero
    assert abs(total_applied + total_reactions) < 1e-3
    
    print("  [PASS] LoadedBeam equilibrium tests passed")


def test_von_mises():
    """Test Von Mises stress calculation."""
    print("Testing von_mises...")
    
    material = Material(name="Steel", E=200e9, G=80e9)
    section = Section(
        name="Rect",
        A=0.01,
        Iy=1e-6,
        Iz=1e-6,
        J=2e-6,
        y_max=0.05,
        z_max=0.05,
    )
    
    supports = [
        Support(x=0.0, type="111111"),
        Support(x=5.0, type="111111"),
    ]
    
    beam = Beam1D(L=5.0, material=material, section=section, supports=supports)
    
    loads = LoadCase(name="Test")
    # Apply combined loading to generate axial, bending, and shear stresses
    loads.add_point_force(PointForce(
        point=np.array([2.5, 0.0, 0.0]),
        force=np.array([1000.0, 1000.0, -1000.0])
    ))
    
    loaded_beam = LoadedBeam(beam, loads)
    
    vm_res = loaded_beam.von_mises(points=11)
    
    assert isinstance(vm_res, Result)
    assert len(list(vm_res)) == 11
    
    # Stress should be non-negative (it's a magnitude)
    assert vm_res.min >= -1e-9 
    
    # Should have some stress
    assert vm_res.max > 0.0
    
    print("  [PASS] von_mises tests passed")


def test_error_cases():
    """Test error handling."""
    print("Testing error cases...")
    
    # Test solve_x_reactions with no supports
    try:
        solve_x_reactions([], LoadCase(name="Test"))
        assert False, "Should have raised ValueError"
    except ValueError:
        pass
    
    # Test solve_transverse_reactions with invalid axis
    material = Material(name="Steel", E=200e9, G=80e9)
    section = Section(name="Rect", A=0.01, Iy=1e-6, Iz=1e-6, J=2e-6, y_max=0.05, z_max=0.05)
    supports = [Support(x=0.0, type="111111")]
    beam = Beam1D(L=5.0, material=material, section=section, supports=supports)
    
    try:
        solve_transverse_reactions(beam, LoadCase(name="Test"), axis="x")
        assert False, "Should have raised ValueError"
    except ValueError:
        pass
    
    print("  [PASS] Error case tests passed")


def run_all_tests():
    """Run all tests."""
    print("=" * 60)
    print("Running analysis.py test suite")
    print("=" * 60)
    print()
    
    tests = [
        test_result_class,
        test_analysis_result,
        test_accumulate_loads,
        test_linear_interpolation,
        test_hermite_displacement,
        test_solve_fem_1d_simple,
        test_solve_x_reactions,
        test_solve_transverse_reactions,
        test_get_all_loads,
        test_loaded_beam_basic,
        test_loaded_beam_equilibrium,
        test_von_mises,
        test_error_cases,
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            test()
            passed += 1
        except Exception as e:
            print(f"  [FAIL] Test failed: {test.__name__}")
            print(f"    Error: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
        print()
    
    print("=" * 60)
    print(f"Test Results: {passed} passed, {failed} failed")
    print("=" * 60)
    
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)

