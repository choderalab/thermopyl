import os
from pkg_resources import resource_filename


def get_fn(name):
    """Get the full path to one of the reference files shipped for testing

    In the source distribution, these files are in ``thermopyl/data``,
    but on installation, they're moved to somewhere in the user's python
    site-packages directory.

    Parameters
    ----------
    name : str
        Name of the file to load (with respect to the reference/ folder).
    """

    fn = resource_filename('thermopyl', os.path.join('data', name))

    if not os.path.exists(fn):
        raise ValueError('Sorry! %s does not exists. If you just '
            'added it, you\'ll have to re install' % fn)

    return fn
