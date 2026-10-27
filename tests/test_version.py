"""The runtime version must match the installed distribution's version.

Wheels are built by hatchling, which does not rewrite versioneer's
`_version.py`; reading that file from an installed wheel yields `0+unknown`.
"""

from importlib.metadata import version

import rtctools_channel_flow


def test_version_matches_distribution():
    assert rtctools_channel_flow.__version__ == version("rtc-tools-channel-flow")
