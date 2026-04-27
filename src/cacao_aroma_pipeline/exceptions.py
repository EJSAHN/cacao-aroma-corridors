class PipelineError(Exception):
    """Base package exception."""


class DetectionError(PipelineError):
    """Raised when a raw file cannot be classified."""


class ParsingError(PipelineError):
    """Raised when a supported file cannot be parsed."""
