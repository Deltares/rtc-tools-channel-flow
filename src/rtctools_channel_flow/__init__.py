# Get version.
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("rtc-tools-channel-flow")
except PackageNotFoundError:
    # Not installed, e.g. imported from a source checkout when building docs.
    from ._version import get_versions

    __version__ = get_versions()["version"]
    del get_versions
