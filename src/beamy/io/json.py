from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict

from beamy.loads.loadcase import LoadCase
from beamy.loads.member_loads import MemberDistributedLoad, MemberPointLoad, MemberPointMoment
from beamy.loads.nodal_loads import NodalLoad
from beamy.loads.supports import PrescribedDisplacement, Support
from beamy.model.frame import Frame
from beamy.model.member import Member
from beamy.model.node import Node

SCHEMA_VERSION = "0.1.0"


def save_frame(frame: Frame, path: str) -> None:
    data = _frame_to_dict(frame)
    _write_json(path, data)


def load_frame(path: str) -> Frame:
    data = _read_json(path)
    return _frame_from_dict(data)


def save_loadcase(loadcase: LoadCase, path: str) -> None:
    data = _loadcase_to_dict(loadcase)
    _write_json(path, data)


def load_loadcase(path: str) -> LoadCase:
    data = _read_json(path)
    return _loadcase_from_dict(data)


def _frame_to_dict(frame: Frame) -> Dict[str, Any]:
    return {
        "schema": SCHEMA_VERSION,
        "nodes": [asdict(node) for node in frame.nodes.values()],
        "members": [
            {
                "id": member.id,
                "end_i": member.end_i,
                "end_j": member.end_j,
                "kind": member.kind,
                "area": member.area,
                "youngs_modulus": member.youngs_modulus,
                "shear_modulus": member.shear_modulus,
                "iy": member.iy,
                "iz": member.iz,
                "j": member.j,
                "density": member.density,
                "stations": list(member.stations),
                "release_i": list(member.release_i),
                "release_j": list(member.release_j),
                "meta": member.meta,
            }
            for member in frame.members.values()
        ],
    }


def _frame_from_dict(data: Dict[str, Any]) -> Frame:
    if "nodes" not in data or "members" not in data:
        raise ValueError("invalid frame json: missing nodes or members")
    frame = Frame()
    for node_data in data["nodes"]:
        meta = node_data["meta"] if "meta" in node_data else {}
        frame.add_node(Node(id=node_data["id"], xyz=node_data["xyz"], meta=meta))
    for member_data in data["members"]:
        density = member_data["density"] if "density" in member_data else 0.0
        shear_modulus = member_data["shear_modulus"] if "shear_modulus" in member_data else 0.0
        iy = member_data["iy"] if "iy" in member_data else 0.0
        iz = member_data["iz"] if "iz" in member_data else 0.0
        j = member_data["j"] if "j" in member_data else 0.0
        stations_raw = member_data["stations"] if "stations" in member_data else []
        release_i_raw = member_data["release_i"] if "release_i" in member_data else [False] * 6
        release_j_raw = member_data["release_j"] if "release_j" in member_data else [False] * 6
        meta_member = member_data["meta"] if "meta" in member_data else {}
        frame.add_member(
            Member(
                id=member_data["id"],
                end_i=member_data["end_i"],
                end_j=member_data["end_j"],
                kind=member_data["kind"],
                area=member_data["area"],
                youngs_modulus=member_data["youngs_modulus"],
                shear_modulus=shear_modulus,
                iy=iy,
                iz=iz,
                j=j,
                density=density,
                stations=tuple(stations_raw),
                release_i=tuple(release_i_raw),  # type: ignore[arg-type]
                release_j=tuple(release_j_raw),  # type: ignore[arg-type]
                meta=meta_member,
            )
        )
    return frame


def _loadcase_to_dict(loadcase: LoadCase) -> Dict[str, Any]:
    return {
        "schema": SCHEMA_VERSION,
        "name": loadcase.name,
        "nodal_loads": [asdict(load) for load in loadcase.nodal_loads],
        "member_loads": [_member_load_to_dict(load) for load in loadcase.member_loads],
        "supports": [asdict(support) for support in loadcase.supports],
        "prescribed_displacements": [asdict(disp) for disp in loadcase.prescribed_displacements],
    }


def _loadcase_from_dict(data: Dict[str, Any]) -> LoadCase:
    if "name" not in data:
        raise ValueError("invalid loadcase json: missing name")
    nodal_loads_raw = data["nodal_loads"] if "nodal_loads" in data else []
    member_loads_raw = data["member_loads"] if "member_loads" in data else []
    supports_raw = data["supports"] if "supports" in data else []
    prescribed_raw = data["prescribed_displacements"] if "prescribed_displacements" in data else []
    nodal_loads = [NodalLoad(**item) for item in nodal_loads_raw]
    member_loads = [_member_load_from_dict(item) for item in member_loads_raw]
    supports = [Support(target=item["target"], restrained_dofs=tuple(item["restrained_dofs"])) for item in supports_raw]  # type: ignore[arg-type]
    prescribed = [PrescribedDisplacement(target=item["target"], values=tuple(item["values"])) for item in prescribed_raw]  # type: ignore[arg-type]
    return LoadCase(
        name=data["name"],
        nodal_loads=nodal_loads,
        member_loads=member_loads,
        supports=supports,
        prescribed_displacements=prescribed,
    )


def _member_load_to_dict(load: Any) -> Dict[str, Any]:
    if isinstance(load, MemberPointLoad):
        return {"type": "point_load", "member_id": load.member_id, "s": load.s, "vector": list(load.vector), "coord_sys": load.coord_sys}
    if isinstance(load, MemberPointMoment):
        return {"type": "point_moment", "member_id": load.member_id, "s": load.s, "vector": list(load.vector), "coord_sys": load.coord_sys}
    if isinstance(load, MemberDistributedLoad):
        return {
            "type": "distributed_load",
            "member_id": load.member_id,
            "s0": load.s0,
            "s1": load.s1,
            "w0": list(load.w0),
            "w1": list(load.w1),
            "coord_sys": load.coord_sys,
        }
    raise ValueError("unknown member load type")


def _member_load_from_dict(item: Dict[str, Any]) -> Any:
    if "type" not in item:
        raise ValueError("member load dict missing type")
    if item["type"] == "point_load":
        return MemberPointLoad(member_id=item["member_id"], s=item["s"], vector=tuple(item["vector"]), coord_sys=item["coord_sys"])  # type: ignore[arg-type]
    if item["type"] == "point_moment":
        return MemberPointMoment(member_id=item["member_id"], s=item["s"], vector=tuple(item["vector"]), coord_sys=item["coord_sys"])  # type: ignore[arg-type]
    if item["type"] == "distributed_load":
        return MemberDistributedLoad(
            member_id=item["member_id"],
            s0=item["s0"],
            s1=item["s1"],
            w0=tuple(item["w0"]),  # type: ignore[arg-type]
            w1=tuple(item["w1"]),  # type: ignore[arg-type]
            coord_sys=item["coord_sys"],
        )
    raise ValueError(f"unknown member load type {item['type']}")


def _write_json(path: str, data: Dict[str, Any]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


def _read_json(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

