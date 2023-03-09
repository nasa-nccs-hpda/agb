
import csv
from pathlib import Path

from osgeo import ogr
from osgeo.osr import SpatialReference

from core.model.Envelope import Envelope


# -----------------------------------------------------------------------------
# class NeonSite
# -----------------------------------------------------------------------------
class NeonSite(object):

    # -------------------------------------------------------------------------
    # __init__
    # -------------------------------------------------------------------------
    def __init__(self):

        self._abbreviation = None
        self._name = None
        self._center = None
        self._envelope = None

    # -------------------------------------------------------------------------
    # abbreviation
    # -------------------------------------------------------------------------
    @property
    def abbreviation(self) -> str:
        return self._abbreviation

    # -------------------------------------------------------------------------
    # abbreviation Setter
    # -------------------------------------------------------------------------
    @abbreviation.setter
    def abbreviation(self, value: str) -> None:
        self._abbreviation = value

    # -------------------------------------------------------------------------
    # name
    # -------------------------------------------------------------------------
    @property
    def name(self) -> str:
        return self._name

    # -------------------------------------------------------------------------
    # name Setter
    # -------------------------------------------------------------------------
    @name.setter
    def name(self, value: str) -> None:
        self._name = value

    # -------------------------------------------------------------------------
    # center
    # -------------------------------------------------------------------------
    @property
    def center(self) -> ogr.Geometry:
        return self._center

    # -------------------------------------------------------------------------
    # setCenter
    # -------------------------------------------------------------------------
    def setCenter(self, utmX: float, utmY: float, epsg: int) -> None:

        self._center = None

        srs = SpatialReference()
        srs.ImportFromEPSG(epsg)

        if self._envelope and \
            not self._envelope.GetSpatialReference().Equal(srs):

            raise RuntimeError('The center and envelope SRSs must match.')

        self._center = ogr.Geometry(ogr.wkbPoint)
        self._center.AssignSpatialReference(srs)
        self._center.AddPoint(utmX, utmY, 0)

    # -------------------------------------------------------------------------
    # envelope
    # -------------------------------------------------------------------------
    @property
    def envelope(self) -> Envelope:
        return self._envelope

    # -------------------------------------------------------------------------
    # setEnvelope
    # -------------------------------------------------------------------------
    def setEnvelope(self, llx: float, lly: float, urx: float, ury: float,
                    epsg: int) -> None:

        self._envelope = None

        srs = SpatialReference()
        srs.ImportFromEPSG(epsg)

        self._envelope = Envelope()
        self._envelope.addPoint(llx, lly, 0.0, srs)
        self._envelope.addPoint(urx, ury, 0.0, srs)

    # -------------------------------------------------------------------------
    # fromCSV
    # -------------------------------------------------------------------------
    def fromCSV(self, siteCsvFileName: Path, siteAbbreviation: str) -> None:

        with open(siteCsvFileName) as csvFile:

            reader = csv.reader(csvFile, quotechar="'")

            for row in reader:
                if row[0] == siteAbbreviation:
                    self._readRow(row)

        if not self._abbreviation:

            raise RuntimeError('Site, ' + \
                               str(siteAbbreviation) + \
                               ', not found in ' + \
                               str(siteCsvFileName))

    # -------------------------------------------------------------------------
    # readRow
    # -------------------------------------------------------------------------
    def _readRow(self, row: str) -> None:

        # Create the spatial reference for the envelope.
        zone = row[7][:-1]
        hemi = row[7][-1].upper()

        # if hemi != 'N' and hemi != 'S':
        #     raise ValueError('Hemisphere' + hemi + ' is neither "N" or "S".')

        if hemi not in ['N', 'S']:
            raise ValueError('Hemisphere' + hemi + ' is neither "N" or "S".')

        geoid = row[4]

        if geoid != 'WGS84':
            raise RuntimeError('Unknown geoid', + geoid)

        self.abbreviation = str(row[0])
        self.name = str(row[1])

        epsg = int('326' + zone)
        self.setCenter(float(row[6]), float(row[5]), epsg)

        lly = float(row[8])
        llx = float(row[9])
        ury = float(row[10])
        urx = float(row[11])
        self.setEnvelope(llx, lly, urx, ury, epsg)

    # -------------------------------------------------------------------------
    # resize
    # -------------------------------------------------------------------------
    def resize(self, lengthInMeters: int) -> None:

        self._envelope = None

        halfLength = int(abs(lengthInMeters) / 2)
        center = self._center.GetPoint()
        lx = center[0] - halfLength
        rx = center[0] + halfLength
        ly = center[1] - halfLength
        uy = center[1] + halfLength
        self._envelope = Envelope()
        self._envelope.addPoint(lx, ly, 0.0, self.center.GetSpatialReference())
        self._envelope.addPoint(rx, uy, 0.0, self.center.GetSpatialReference())
