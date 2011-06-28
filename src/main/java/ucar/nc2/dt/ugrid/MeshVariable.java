/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package ucar.nc2.dt.ugrid;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import ucar.ma2.Array;
import ucar.ma2.DataType;
import ucar.ma2.IndexIterator;
import ucar.ma2.InvalidRangeException;
import ucar.ma2.MAMath;
import ucar.ma2.MAMath.MinMax;
import ucar.ma2.Range;
import ucar.nc2.Attribute;
import ucar.nc2.Dimension;
import ucar.nc2.Variable;
import ucar.nc2.dataset.CoordinateAxis;
import ucar.nc2.dataset.CoordinateSystem;
import ucar.nc2.dataset.NetcdfDataset;
import ucar.nc2.dataset.VariableDS;
import ucar.nc2.dataset.VariableEnhanced;
import ucar.nc2.dt.UGridDataset.Meshset;
import ucar.nc2.dt.UGridDatatype;
import ucar.nc2.dt.ugrid.geom.LatLonPoint2D;
import ucar.nc2.dt.ugrid.geom.LatLonPolygon2D;
import ucar.nc2.dt.ugrid.geom.LatLonRectangle2D;
import ucar.nc2.dt.ugrid.utils.NcdsFactory;
import ucar.nc2.dt.ugrid.utils.NcdsFactory.NcdsTemplate;
import ucar.unidata.geoloc.LatLonPoint;
import ucar.unidata.geoloc.LatLonRect;
import ucar.unidata.geoloc.ProjectionImpl;
import ucar.unidata.util.Format;

/**
 *
 * @author Kyle
 */
public class MeshVariable implements UGridDatatype {

  private VariableDS vs;
  private UGridDataset dataset;
  private Meshset meshset;
  private List<Dimension> mydims;
  private String cellLocation;

  public MeshVariable(UGridDataset dataset, VariableDS vs, Meshset meshset) {
    this.vs = vs;
    this.meshset = meshset;
    this.dataset = dataset;
    this.mydims = vs.getDimensions();
    this.cellLocation = vs.findAttributeIgnoreCase("location").getStringValue();
  }

  public String getName() {
    return vs.getName();
  }

  public VariableDS getConnectivityVariable() {
    return (VariableDS)meshset.getMesh().getConnectivityVariables().get(0);
  }

  public String getNameEscaped() {
    return vs.getNameEscaped();
  }

  public String getDescription() {
    return vs.getDescription();
  }

  public String getUnitsString() {
    return vs.getUnitsString();
  }

  public DataType getDataType() {
    return vs.getDataType();
  }

  public int getRank() {
    return vs.getRank();
  }

  public int[] getShape() {
    int[] shape = new int[mydims.size()];
    for (int i = 0; i < mydims.size(); i++) {
      Dimension d = mydims.get(i);
      shape[i] = d.getLength();
    }
    return shape;
  }

  public List<Attribute> getAttributes() {
    return vs.getAttributes();
  }

  public Attribute findAttributeIgnoreCase(String name) {
    return vs.findAttributeIgnoreCase(name);    
  }

  public String findAttValueIgnoreCase(String attName, String defaultValue) {
    return dataset.getNetcdfDataset().findAttValueIgnoreCase((Variable) vs, attName, defaultValue);
  }

  public List<Dimension> getDimensions() {
    return mydims;
  }

  public Dimension getDimension(int i) {
    return mydims.get(i);
  }

  public Dimension getTimeDimension() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  public Dimension getZDimension() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  public Dimension getYDimension() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  public Dimension getXDimension() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  public Dimension getEnsembleDimension() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  public Dimension getRunTimeDimension() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  public int getTimeDimensionIndex() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  public int getZDimensionIndex() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  public int getYDimensionIndex() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  public int getXDimensionIndex() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  public int getEnsembleDimensionIndex() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  public int getRunTimeDimensionIndex() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  public Meshset getMeshset() {
    return meshset;
  }

  public ProjectionImpl getProjection() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  public boolean hasMissingData() {
    return vs.hasMissing();
  }

  public boolean isMissingData(double val) {
    return vs.isMissing(val);
  }

  public MinMax getMinMaxSkipMissingData(Array a) {
    if (!hasMissingData())
      return MAMath.getMinMax(a);

    IndexIterator iter = a.getIndexIterator();
    double max = -Double.MAX_VALUE;
    double min = Double.MAX_VALUE;
    while (iter.hasNext()) {
      double val = iter.getDoubleNext();
      if (isMissingData(val))
        continue;
      if (val > max)
        max = val;
      if (val < min)
        min = val;
    }
    return new MAMath.MinMax(min, max);
  }

  public float[] setMissingToNaN(float[] values) {
    if (!vs.hasMissing()) return values;
    final int length = values.length;
    for (int i = 0; i < length; i++) {
      double value = values[i];
      if (vs.isMissing(value))
        values[i] = Float.NaN;
    }
    return values;
  }

  public Array readDataSlice(int rt_index, int e_index, int t_index, int z_index, int y_index, int x_index) throws IOException {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  public Array readDataSlice(int t_index, int z_index, int y_index, int x_index) throws IOException {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  public Array readVolumeData(int t_index) throws IOException {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  public UGridDatatype makeSubset(Range rt_range, Range e_range, Range t_range, Range z_range, Range y_range, Range x_range) throws InvalidRangeException {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  public UGridDatatype makeSubset(Range t_range, Range z_range, LatLonRect bbox, int z_stride, int y_stride, int x_stride) throws InvalidRangeException {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  public UGridDataset subset(LatLonRect bounds) {
    LatLonRectangle2D r = new LatLonRectangle2D(new LatLonPoint2D.Double(bounds.getUpperLeftPoint().getLatitude(), bounds.getUpperLeftPoint().getLongitude()), new LatLonPoint2D.Double(bounds.getLowerRightPoint().getLatitude(), bounds.getLowerRightPoint().getLongitude()));
    LatLonPolygon2D p = new LatLonPolygon2D.Double(r);
    List<Cell> containedCells = meshset.getMesh().getCellsInPolygon(p);

    // Create a new subsat UGridDataset and return
    try {
      NetcdfDataset ncd = NcdsFactory.getNcdsFromTemplate(NcdsTemplate.UGRID);

      for (Attribute a : dataset.getGlobalAttributes()){
        ncd.addAttribute(null, a);
      }

      for (CoordinateSystem cs : dataset.getNetcdfDataset().getCoordinateSystems()) {
        for (Dimension d : cs.getDomain()) {
          if (ncd.findDimension(d.getName()) == null) {
            ncd.addDimension(null, d);
            Variable vd = dataset.getNetcdfFile().findVariable(d.getName());
            if (vd != null) {
              ncd.addVariable(null, new VariableDS(null, vd, true));
            }
          }
          ncd.finish();
        }
        ncd.addCoordinateSystem(cs);
        ncd.finish();
      }

      for (CoordinateAxis ax : dataset.getNetcdfDataset().getCoordinateAxes()) {
        ncd.addCoordinateAxis(ax);
        ncd.finish();
      }

      // Variable describing the Mesh (Mesh1, Mesh2, etc)
      ncd.addVariable(null, meshset.getDescriptionVariable());

      // Connectivity Variable for this MeshVariable
      ncd.addVariable(null, getConnectivityVariable());
      
      ncd.addVariable(null, vs);

      // Now add the data
      //vs.setCachedData(null);
      
      ncd.finish();
      
      return new UGridDataset(ncd);
    } catch (URISyntaxException e) {
      System.out.println(e);
    } catch (FileNotFoundException e) {
      System.out.println(e);
    } catch (IOException e) {
      System.out.println(e);
    }
    return null;
  }

  public double readPointData(LatLonPoint point) throws IOException {
    // Find the closest R-Tree Cell
    final LatLonPoint2D p = new LatLonPoint2D.Double(point.getLatitude(), point.getLongitude());
    Cell c = meshset.getMesh().getCellFromLatLon(p);
    double z = -1;

    List<? extends Entity> e;
    if (cellLocation.equals("node")) {
      e = c.getNodes();
    } else if (cellLocation.equals("face")) {
      e = c.getFaces();
    } else if (cellLocation.equals("edge")) {
      e = c.getEdges();
    } else {
      e = null;
    }
    if ((e != null) && (e.size() > 0)) {

      // Sort the collection of Entities by distance from the query point.
      // This should offer different ways to calculate the distance of
      // the closest point.

      // "e" is all of the Entities in the Cell that are on the variable's
      // location (Node OR Edge OR Face).

      Collections.sort(e, new Comparator() {

        public int compare(Object o1, Object o2) {
          Entity e1 = (Entity) o1;
          Entity e2 = (Entity) o2;
          if (e1.getGeoPoint().distance(p) == e2.getGeoPoint().distance(p)) {
            return 0;
          } else if (e1.getGeoPoint().distance(p) > e2.getGeoPoint().distance(p)) {
            return 1;
          } else {
            return -1;
          }
        }
      });

      // Get the closest Entities DataIndex into the NetCDF file.
      int in = e.get(0).getDataIndex();
      try {

        // Need to compute actual ranges here, not assume it is (time,z,entity)
        List<Range> r = new ArrayList<Range>();
        // Time (first)
        r.add(new Range(0, 0));
        // Sigma (first)
        r.add(new Range(0, 0));
        // Data (DataIndex from Cell)
        r.add(new Range(in, in));

        float[] ret1D = (float[]) vs.read(r).copyTo1DJavaArray();
        z = ret1D[0];
      } catch (InvalidRangeException ex) {
        System.out.println(ex);
      }
    }
    return z;
  }

  public String getInfo() {
    StringBuilder buf = new StringBuilder(200);
    buf.setLength(0);
    buf.append(getName());
    Format.tab(buf, 30, true);
    buf.append(getUnitsString());
    Format.tab(buf, 60, true);
    buf.append(hasMissingData());
    Format.tab(buf, 66, true);
    buf.append(getDescription());
    return buf.toString();
  }

  public VariableDS getVariable() {
    return vs;
  }

  public int compareTo(UGridDatatype g) {
    return getName().compareTo(g.getName());
  }

}