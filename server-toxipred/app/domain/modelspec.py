from dataclasses import dataclass
from typing import Literal, Optional
from pathlib import Path

ModelKind = Literal["pickle"]

@dataclass(frozen=True)
class ModelSpec:
    path: Path
    kind: ModelKind
    features: Optional[list[str]] = None
    note: str = ""