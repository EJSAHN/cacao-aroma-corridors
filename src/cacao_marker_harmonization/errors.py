class DataError(ValueError):
    """An input or configuration fails an explicit analysis requirement."""

class HeaderNotFound(DataError):
    """A required worksheet header is absent."""
