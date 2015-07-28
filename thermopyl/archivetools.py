import os
import tarfile
import feedparser
from .utils import make_path
from six.moves import urllib_parse, urllib

THERMOML_FEEDS = {
"jced":"http://trc.nist.gov/RSS/jced.xml",
"jct":"http://trc.nist.gov/RSS/jct.xml",
"fpe":"http://trc.nist.gov/RSS/fpe.xml",
"tca":"http://trc.nist.gov/RSS/tca.xml",
"ijt":"http://trc.nist.gov/RSS/ijt.xml"
}

THERMOML_TARBALLS = {
"jced": "http://trc.nist.gov/ThermoML/JCED.tgz",
"jct": "http://trc.nist.gov/ThermoML/JCT.tgz",
"fpe": "http://trc.nist.gov/ThermoML/FPE.tgz",
"tca": "http://trc.nist.gov/ThermoML/TCA.tgz",
"ijt": "http://trc.nist.gov/ThermoML/IJT.tgz"
}


def update_archive(thermoml_path=None):
    """Use RSS feeds to find and download ThermoML tar files
    from the ThermoML archive, then download any missing entries by enumerating the
    RSS feeds.  The output will be a flat directory of XML files in `thermoml_path`

    Parameters
    ----------
    thermoml_path : str, optional, default=None
        If specified, use this path to store ThermoML XML files.  If None,
        use the THERMOML_PATH environment variable.
    """
    if thermoml_path is None:
        # Try to obtain the path to the local ThermoML Archive mirror from an environment variable.
        try:
            # Check THERMOML_PATH environment variable
            thermoml_path = os.environ["THERMOML_PATH"]
        except:
            # Use default path of ~/.thermoml
            thermoml_path = os.path.join(os.environ["HOME"], '.thermoml')

    if not os.path.exists(thermoml_path):
        os.makedirs(thermoml_path)

    for key, url in THERMOML_TARBALLS.items():
        print("Downloading %s" % url)
        file = urllib.request.URLopener()
        local_filename = key + ".tgz"
        file.retrieve(url, local_filename)
        tarball = tarfile.open(local_filename)
        tarball.extractall(thermoml_path)

    # Update local repository according to feeds.
    for key, url in THERMOML_FEEDS.items():
        print("Fetching RSS %s" % url)
        feed = feedparser.parse(url)
        for entry in feed["entries"]:
            link = entry["link"]
            base_filename = urllib_parse.urlsplit(link).path
            base_filename = os.path.split(base_filename)[-1]  # Flattens the directory structure to a flat directory to make the tar and RSS files compatible.
            filename = os.path.join(thermoml_path, base_filename)
            make_path(filename)
            if os.path.exists(filename):
                print("Already downloaded %s from %s" % (filename, link))
            else:
                print("Fetching %s from %s" % (filename, link))
                urllib.request.urlretrieve (link, filename)
