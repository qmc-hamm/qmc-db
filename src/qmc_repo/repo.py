"""Repository access for QMC configuration datasets."""

import intake
import s3fs


class Repo:
    """Interface to the QMC configurations repository stored on S3.

    Provides access to an intake catalog of QMC configuration data files
    hosted on the Open Storage Network (OSN). The catalog entries can be
    accessed as attributes on this class.

    Attributes:
        endpoint_url: The S3 endpoint URL for OSN.
        catalog: The intake catalog containing dataset entries.
        fs: An s3fs filesystem for accessing HDF5 files directly.
    """

    def __init__(self):
        self.endpoint_url = "https://uri.osn.mghpcc.org"

        self.catalog = intake.open_catalog(
            's3://phy240060/qmchamm.yaml',
            storage_options={
                'anon': True,
                'endpoint_url': self.endpoint_url}
            )

        # Initialize s3fs to allow access to hdf5 files in S3
        self.fs = s3fs.S3FileSystem(
            endpoint_url=self.endpoint_url,
            anon=True
        )

    @property
    def entries(self):
        """List of available dataset entry names in the catalog."""
        return [c for c in self.catalog]

    def get_source(self, name):
        """Get a catalog entry by name.

        Args:
            name: The name of the catalog entry to retrieve.

        Returns:
            The intake catalog source for the given entry.
        """
        return getattr(self.catalog, name)

    def __getattr__(self, name):
        """Delegate attribute access to the underlying catalog.

        Allows catalog entries to be accessed directly as attributes on
        the Repo instance (e.g., repo.hydrogen instead of repo.catalog.hydrogen).

        Args:
            name: The attribute name to look up.

        Returns:
            The corresponding attribute from the catalog.
        """
        return getattr(self.catalog, name)
