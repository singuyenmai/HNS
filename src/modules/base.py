from abc import ABC, abstractmethod
from typing import Any

class PipelineStep(ABC):
    """
    Abstract base class for all modules in the H-NS extraction pipeline.
    """
    @abstractmethod
    def run(self, *args: Any, **kwargs: Any) -> Any:
        """
        Execute the pipeline step.
        """
        pass
