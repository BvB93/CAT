"""Test the :mod:`sphinx` documentation generation."""

from os.path import join
from typing import Optional

from sphinx.application import Sphinx
from nanoutils import delete_finally, ignore_if

try:
    import nanoCAT  # noqa: F401
    NANOCAT_EX: Optional[ImportError] = None
except ImportError as ex:
    NANOCAT_EX = ex

SRCDIR = CONFDIR = 'docs'
OUTDIR = join('tests', 'test_files', 'build')
DOCTREEDIR = join('tests', 'test_files', 'build', 'doctrees')


@ignore_if(NANOCAT_EX)
@delete_finally(OUTDIR)
def test_sphinx_build() -> None:
    """Test :meth:`sphinx.application.Sphinx.build`."""
    app = Sphinx(SRCDIR, CONFDIR, OUTDIR, DOCTREEDIR, buildername='html', warningiserror=True)
    app.build(force_all=True)
