# Note that we need to fall back to the hard-coded version if either
# setuptools_scm can't be imported or setuptools_scm can't determine the
# version, so we catch the generic 'Exception'.
try:
    from setuptools_scm import get_version
    version = get_version(root='..', relative_to=__file__)
except Exception:
<<<<<<< HEAD
    version = '0.1.dev31+g4f55b42'
=======
    version = '0.3.dev13+g1061680.d20230308'
>>>>>>> 8bbb7e41ca4c73db17f3186a2a042d38198c0d45
