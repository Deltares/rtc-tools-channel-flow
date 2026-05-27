This is the ChannelFlow library for RTC-Tools 2.0.  It contains both simple as well as hydraulic routing models.

Visit the RTC-Tools website at:

    https://www.deltares.nl/en/software/rtc-tools/

The documentation of the RTC-Tools channel flow library can be found at 
	
	https://rtc-tools-channel-flow.readthedocs.io/
	
Note that this documentation is still under development and not yet complete. 	
    
## Documentation

Documentation and examples can be found on [readthedocs](https://rtc-tools-simulation.readthedocs.io).

To build the documentation, the required dependencies are in the `docs` dependency group from `pyproject.toml`.
Run these commands from the repository root:

```bash
uv sync --group docs
uv run sphinx-build -b html docs/source docs/build/html
```

The built HTML pages will be in `docs/build/html/`.