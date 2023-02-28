import os
import unittest

from osgeo import ogr
from osgeo.osr import SpatialReference

from agb.model.NeonSite import NeonSite
from core.model.Envelope import Envelope


# -----------------------------------------------------------------------------
# class NeonSiteTestCase
#
# python -m unittest discover agb/model/tests/
# python -m unittest agb.model.tests.test_NeonSite
# python -m unittest agb.model.tests.test_NeonSite.NeonSiteTestCase.testCenter
# -----------------------------------------------------------------------------
class NeonSiteTestCase(unittest.TestCase):

    # -------------------------------------------------------------------------
    # testReadRow
    # -------------------------------------------------------------------------
    def testReadRow(self):

        row = ['MLBS',
               'Mountain Lake Biological Station NEON',
           	   '37.37831',
               '-80.5248',
               'WGS84',
               '4136943',
               '542067.6',
               '17N',
               '4136443',
               '541567.6',
               '4137443',
               '542567.6']

        site = NeonSite()
        site._readRow(row)

    # -------------------------------------------------------------------------
    # testInit
    # -------------------------------------------------------------------------
    def testInit(self):

        site = NeonSite()
        self.assertIsNotNone(site)

    # -------------------------------------------------------------------------
    # testAbbreviation
    # -------------------------------------------------------------------------
    def testAbbreviation(self):

        site = NeonSite()
        site.abbreviation = 'ABCD'
        self.assertEqual(site.abbreviation, 'ABCD')

    # -------------------------------------------------------------------------
    # testCenter
    # -------------------------------------------------------------------------
    def testCenter(self):

        site = NeonSite()
        site.setCenter(542067.6, 4136943, 32617)
        self.assertIsNotNone(site.center)
        self.assertIsInstance(site.center, ogr.Geometry)

    # -------------------------------------------------------------------------
    # testEnvelope
    # -------------------------------------------------------------------------
    def testEnvelope(self):

        site = NeonSite()
        site.setCenter(542067.6, 4136943, 32617)

        lly = 4136443
        llx = 541567.6
        ury = 4137443
        urx = 542567.6
        site.setEnvelope(llx, lly, urx, ury, 32617)

        self.assertIsNotNone(site.envelope)
        self.assertEqual(site.envelope.ulx(), llx)
        self.assertEqual(site.envelope.uly(), ury)
        self.assertEqual(site.envelope.lrx(), urx)
        self.assertEqual(site.envelope.lry(), lly)

    # -------------------------------------------------------------------------
    # testFromCSV
    # -------------------------------------------------------------------------
    def testFromCSV(self):

        csvName = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             'NeonSites.csv')

        abby = NeonSite()
        abby.fromCSV(csvName, 'ABBY')
        self.assertEqual(abby.name, 'Abby Road NEON')

        yell = NeonSite()
        yell.fromCSV(csvName, 'YELL')
        self.assertEqual(yell.name, 'Yellowstone National Park NEON')

        srs = SpatialReference()
        srs.ImportFromEPSG(32612)
        env = Envelope()
        env.addPoint(535852.2, 4977386, 0.0, srs)
        env.addPoint(536852.2, 4978386, 0.0, srs)
        self.assertEqual(yell.envelope.lrx(), env.lrx())
        self.assertEqual(yell.envelope.lry(), env.lry())
        self.assertEqual(yell.envelope.ulx(), env.ulx())
        self.assertEqual(yell.envelope.uly(), env.uly())

        centerX = yell.envelope.ulx() + \
                  (yell.envelope.lrx() - yell.envelope.ulx()) / 2

        centerY = yell.envelope.lry() + \
                  (yell.envelope.uly() - yell.envelope.lry()) / 2

        center = ogr.Geometry(ogr.wkbPoint)
        center.AssignSpatialReference(srs)
        center.AddPoint(centerX, centerY, 0)
        self.assertTrue(yell.center.Equal(center))

    # -------------------------------------------------------------------------
    # testResize
    # -------------------------------------------------------------------------
    def testResize(self):

        csvName = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               'NeonSites.csv')

        abby = NeonSite()
        abby.fromCSV(csvName, 'ABBY')

        resizeAmount = 10
        centerX = 552075.5
        centerY = 5067870.0
        expUlx = centerX - resizeAmount / 2
        expUly = centerY + resizeAmount / 2
        expLrx = centerX + resizeAmount / 2
        expLry = centerY - resizeAmount / 2
        abby.resize(resizeAmount)
        self.assertEqual(abby.envelope.ulx(), expUlx)
        self.assertEqual(abby.envelope.uly(), expUly)
        self.assertEqual(abby.envelope.lrx(), expLrx)
        self.assertEqual(abby.envelope.lry(), expLry)
