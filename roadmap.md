# ROADMAP

## Needed for 4.0
### Still to do
- Up to date doc (for install and use)
- No documentation build by default
- code beautifying (linter)
- Better test for TXT
- Stamped TXT set
- Merging Boris changes
- Merging Natasha changes
- Checksum the distributed particle set

### Maybe
- THERMUS and THERMUS_LIB
Now we have two env variables: THERMUS and THERMUS_LIB
As we can't predict if libs will be in /lib or /lib64.
We could either force install in /lib, but that would be malpracice on Linux
or keep those two variables.
We could also hardcode THERMUS_LIB in use_thermus.C, making it reusable on a given platform, for a different THERMUS path.
For now, kept as it is.

### Already done:
- Fix rdict.pcm install
- Fix particle install
- LaTeX doc working (on Linux)
- removed rootlogon (replaced by use_thermus.C)
- Run on MacOS (without LaTeX doc generation)
- Autotests by githlab workflow
- Library path unders Thermus
- Scripts in their own subdir
- Removed support for Makefile build


## Before next upstream push
### Already done
- build on macos (without doc)
- roadmap as an md file


## Later
- Fix xxHash test
- TXT files optimisation (Maybe? Not now)
- Local TXT files?
- Automatic PDG upgrade?




