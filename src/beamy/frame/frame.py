from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Dict, Set, Tuple, Optional, TYPE_CHECKING
import numpy as np
from .node import Node
from .member import Member

if TYPE_CHECKING:
    from ..core.loads import LoadCase
    from .analysis import FrameAnalysisResult, FrameAnalysisSettings
    from .analysis import _FrameSolveState
    from .results import MemberDemandProvider, MemberDemand

@dataclass
class Frame:
    """A collection of nodes and members forming a 3D structural system."""
    members: List[Member]
    nodes: Dict[str, Node] = field(init=False, default_factory=dict)

    # Cached analysis state (populated by analyze(); does not mutate geometry)
    _analysis_result: Optional["FrameAnalysisResult"] = field(init=False, default=None, repr=False)
    _demand_provider: Optional["MemberDemandProvider"] = field(init=False, default=None, repr=False)
    _solve_state: Optional["_FrameSolveState"] = field(init=False, default=None, repr=False)
    _last_load_case: Optional["LoadCase"] = field(init=False, default=None, repr=False)
    _last_settings: Optional["FrameAnalysisSettings"] = field(init=False, default=None, repr=False)
    
    def __post_init__(self) -> None:
        raise NotImplementedError("Use Frame.from_members(members) or Frame.from_nodes_and_members(nodes, members)")
    
    @classmethod
    def from_members(cls, members: List[Member], node_tolerance: float = 1e-6) -> Frame:
        """Build a Frame from members, auto-generating nodes at member endpoints.
        
        Coincident endpoints (within node_tolerance) are merged into single nodes.
        Node IDs are auto-generated as N0, N1, N2, ...
        
        Args:
            members: List of Member objects (with start/end positions)
            node_tolerance: Distance threshold for merging coincident nodes
        
        Returns:
            Frame with auto-generated nodes and connectivity
        """
        frame = cls.__new__(cls)
        frame.members = members
        frame.nodes = {}
        frame._analysis_result = None
        frame._demand_provider = None
        frame._solve_state = None
        frame._last_load_case = None
        frame._last_settings = None
        
        # Collect all unique positions and build nodes
        positions: List[np.ndarray] = []
        position_to_node_id: Dict[int, str] = {}
        
        for member in members:
            for pos in [member.start, member.end]:
                # Check if this position is close to an existing one
                found_idx = None
                for idx, existing_pos in enumerate(positions):
                    if np.linalg.norm(pos - existing_pos) < node_tolerance:
                        found_idx = idx
                        break
                
                if found_idx is None:
                    # New unique position
                    node_id = f"N{len(positions)}"
                    positions.append(pos.copy())
                    position_to_node_id[len(positions) - 1] = node_id
                    frame.nodes[node_id] = Node(id=node_id, position=pos.copy())
        
        # Now assign node IDs to members by finding closest nodes
        for member in members:
            # Find node for start
            start_idx = None
            for idx, pos in enumerate(positions):
                if np.linalg.norm(member.start - pos) < node_tolerance:
                    start_idx = idx
                    break
            if start_idx is None:
                raise RuntimeError(f"Member {member.id}: could not find node for start position")
            
            # Find node for end
            end_idx = None
            for idx, pos in enumerate(positions):
                if np.linalg.norm(member.end - pos) < node_tolerance:
                    end_idx = idx
                    break
            if end_idx is None:
                raise RuntimeError(f"Member {member.id}: could not find node for end position")
            
            # Set the internal node references
            start_node = frame.nodes[position_to_node_id[start_idx]]
            end_node = frame.nodes[position_to_node_id[end_idx]]
            member.set_nodes(start_node, end_node)
            
            # Update connected_members
            if member.id not in start_node.connected_members:
                start_node.connected_members.append(member.id)
            if member.id not in end_node.connected_members:
                end_node.connected_members.append(member.id)
        
        # Validate
        frame._validate_member_ids()
        frame._validate_supports()
        return frame
    
    @classmethod
    def from_nodes_and_members(cls, nodes: List[Node], members: List[Member]) -> Frame:
        """Build a Frame from explicit nodes and members.
        
        This is the legacy low-level constructor. Members must reference node positions
        that match the provided nodes (within tolerance).
        
        Args:
            nodes: List of Node objects with explicit IDs
            members: List of Member objects (with start/end positions)
        
        Returns:
            Frame with explicit node connectivity
        """
        frame = cls.__new__(cls)
        frame.members = members
        frame.nodes = {}
        frame._analysis_result = None
        frame._demand_provider = None
        frame._solve_state = None
        frame._last_load_case = None
        frame._last_settings = None
        
        # Build node dict
        for n in nodes:
            if n.id in frame.nodes:
                raise ValueError(f"Duplicate node ID: {n.id}")
            frame.nodes[n.id] = n
        
        # Match members to nodes by position
        node_tolerance = 1e-6
        for member in members:
            # Find matching start node
            start_node = None
            for node in frame.nodes.values():
                if np.linalg.norm(member.start - node.position) < node_tolerance:
                    start_node = node
                    break
            if start_node is None:
                raise ValueError(f"Member {member.id}: no node found at start position {member.start}")
            
            # Find matching end node
            end_node = None
            for node in frame.nodes.values():
                if np.linalg.norm(member.end - node.position) < node_tolerance:
                    end_node = node
                    break
            if end_node is None:
                raise ValueError(f"Member {member.id}: no node found at end position {member.end}")
            
            # Set nodes and update connectivity
            member.set_nodes(start_node, end_node)
            if member.id not in start_node.connected_members:
                start_node.connected_members.append(member.id)
            if member.id not in end_node.connected_members:
                end_node.connected_members.append(member.id)
        
        frame._validate_member_ids()
        frame._validate_supports()
        return frame
    
    def _validate_member_ids(self) -> None:
        """Validate that member IDs are unique and members have valid length."""
        mids = set()
        for m in self.members:
            if m.id in mids:
                raise ValueError(f"Duplicate member ID: {m.id}")
            mids.add(m.id)
            if m.length < 1e-10:
                raise ValueError(f"Member {m.id}: zero-length")
    
    def _validate_supports(self) -> None:
        """Validate that the frame has sufficient support constraints for stability."""
        constrained_dofs: Set[Tuple[str, int]] = set()

        # Count node supports
        for node in self.nodes.values():
            if node.support:
                for d in range(6):
                    if node.support[d] == "1":
                        constrained_dofs.add((node.id, d))

        # Count member end constraints
        for m in self.members:
            if m.constraints:
                start_node_id = m._start_node.id if m._start_node else None
                end_node_id = m._end_node.id if m._end_node else None
                if start_node_id:
                    for d in range(6):
                        if m.constraints[d] == "1":
                            constrained_dofs.add((start_node_id, d))
                if end_node_id:
                    for d in range(6):
                        if m.constraints[6 + d] == "1":
                            constrained_dofs.add((end_node_id, d))

        if not constrained_dofs:
            raise ValueError("No supports found (neither Node.support nor Member.constraints)")
        if len(constrained_dofs) < 6:
            raise ValueError("Insufficient supports for 3D stability (< 6 constrained DOFs)")
    
    @property
    def node_positions(self) -> Dict[str, np.ndarray]: return {nid: n.position for nid, n in self.nodes.items()}
    @property
    def supported_nodes(self) -> List[Node]: return [n for n in self.nodes.values() if n.support]
    @property
    def member_lengths(self) -> Dict[str, float]: return {m.id: m.length for m in self.members}
    
    def get_member(self, mid: str) -> Member:
        for m in self.members:
            if m.id == mid: return m
        raise KeyError(f"Member '{mid}' not found")
    
    def get_node(self, nid: str) -> Node:
        if nid not in self.nodes: raise KeyError(f"Node '{nid}' not found")
        return self.nodes[nid]
    
    def members_at_node(self, nid: str) -> List[Member]:
        return [self.get_member(mid) for mid in self.get_node(nid).connected_members]
    
    def __repr__(self) -> str: return f"Frame({len(self.nodes)} nodes, {len(self.members)} members)"

    def analyze(self, load_case: "LoadCase", settings: Optional["FrameAnalysisSettings"] = None) -> "FrameAnalysisResult":
        """Analyze this frame under a single load case and cache results on the frame.

        Notes:
            - This does NOT mutate the user-defined geometry/topology (nodes/members).
            - The solver may internally auto-split members at attachment points (loads/supports),
              but that expanded model remains internal.
        """
        from .analysis import FrameAnalysisSettings, _solve_frame_internal  # local import avoids circular deps

        use_settings = settings if settings is not None else FrameAnalysisSettings()
        solve_state = _solve_frame_internal(self, load_case=load_case, settings=use_settings)

        # Cache user-facing state
        self._solve_state = solve_state
        self._analysis_result = solve_state.result
        self._demand_provider = solve_state.result.demand_provider
        self._last_load_case = load_case
        self._last_settings = use_settings

        return solve_state.result

    @property
    def analysis_result(self) -> "FrameAnalysisResult":
        if self._analysis_result is None:
            raise RuntimeError("Frame has not been analyzed yet. Call frame.analyze(load_case) first.")
        return self._analysis_result

    @property
    def demand_provider(self) -> "MemberDemandProvider":
        if self._demand_provider is None:
            raise RuntimeError("Demand provider is available after frame.analyze(load_case) completes.")
        return self._demand_provider

    def member_demand(self, member_id: str) -> "MemberDemand":
        from .results import MemberDemand

        return MemberDemand(member_id=member_id, _provider=self.demand_provider)

    @property
    def nodal_displacements(self) -> Dict[str, np.ndarray]:
        """Solved nodal displacement vectors [UX, UY, UZ, RX, RY, RZ] keyed by node_id."""
        return self.analysis_result.nodal_displacements

    @property
    def reactions(self) -> Dict[str, np.ndarray]:
        """Solved support reactions [FX, FY, FZ, MX, MY, MZ] keyed by node_id."""
        return self.analysis_result.reactions

    @property
    def member_end_forces(self) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
        """Solved local member end forces keyed by member_id: (start_forces, end_forces)."""
        return self.analysis_result.member_end_forces
