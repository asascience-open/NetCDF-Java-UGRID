/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ucar.nc2.dt.ugrid;

import cern.colt.list.IntArrayList;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Properties;
import java.util.Set;
import ucar.nc2.Attribute;
import ucar.nc2.Variable;
import ucar.nc2.dataset.CoordinateSystem;
import ucar.nc2.dataset.NetcdfDataset;
import ucar.nc2.dataset.VariableEnhanced;
import ucar.nc2.dt.ugrid.geom.LatLonPoint2D;
import ucar.nc2.dt.ugrid.geom.LatLonPolygon2D;
import ucar.nc2.dt.ugrid.geom.LatLonRectangle2D;
import ucar.nc2.dt.ugrid.rtree.RTree;
import ucar.unidata.geoloc.LatLonPoint;
import ucar.unidata.geoloc.LatLonPointImpl;
import ucar.unidata.geoloc.LatLonRect;

/**
 *
 * @author Kyle
 */
public class Mesh {

  private String name;
  private RTree rtree;
  private List<Cell> cells = new ArrayList<Cell>();
  // A Mesh should only have one connectivity array!
  private ConnectivityVariable connectivity_variable;
  private List<String> locations;
  private List<CoordinateSystem> coordinate_systems = new ArrayList<CoordinateSystem>();
  
  // Standards
  private final static String LOCATIONS_ATTRIBUTE = "locations";
  private final static String COORDINATES = "coordinates";
  private final static String CONNECTIVITY = "connectivity";

  public Mesh(NetcdfDataset ds, VariableEnhanced v) {
    name = v.getName();

    Properties props = new Properties();
    props.setProperty("MaxNodeEntries", "30");
    props.setProperty("MinNodeEntries", "15");
    rtree = new RTree(props);

    processConnectivityVariables(ds, v);
  }

  private void processConnectivityVariables(NetcdfDataset ncd, VariableEnhanced v) {
    ArrayList<Attribute> foundCoords;
    Attribute att = null;
    if (v.findAttributeIgnoreCase(LOCATIONS_ATTRIBUTE) != null) {
      String[] refs = v.findAttributeIgnoreCase(LOCATIONS_ATTRIBUTE).getStringValue().split(" ");
      foundCoords = new ArrayList<Attribute>(refs.length);
      locations = new ArrayList<String>(refs.length);
      Attribute conn = null;
      Variable connection_variable;
      connectivity_variable = null;
      for (String prefix : refs) {
        att = v.findAttributeIgnoreCase(prefix + "_" + COORDINATES);
        conn = v.findAttributeIgnoreCase(prefix + "_" + CONNECTIVITY);
        if (att != null && conn != null && !att.getStringValue().isEmpty() && !conn.getStringValue().isEmpty()) {
          connection_variable = ncd.findVariable(conn.getStringValue());
          // Add support for "edge" cases where connectivity gives the edges,
          // and there is another map to the nodes.
          if (prefix.equalsIgnoreCase("node") && connectivity_variable == null) {
            connectivity_variable = new ConnectivityVariable(connection_variable);
          }
          locations.add(prefix);
          foundCoords.add(att);
        } else {
          System.out.println("Could not find coordinate and connectivity information for this Mesh.");
        }
      }

      String attString;
      findCoord : for (Attribute attr : foundCoords) {
        attString = attr.getStringValue().toLowerCase();
        for (CoordinateSystem coord : ncd.getCoordinateSystems()) {
          if (attString.contains(coord.getLatAxis().getShortName()) &&
                  attString.contains(coord.getLonAxis().getShortName())) {
            coordinate_systems.add(coord);
            break;
          }
        }
      }
      cells = getConnectivityVariable().createCells(locations, coordinate_systems);
    } else {
      System.out.println("No 'locations' attribute on the Mesh.");
    }
  }

  public void buildRTree() {
    for (int i = 0; i < cells.size(); i++) {
      rtree.add(cells.get(i).getPolygon(), i);
    }
  }

  public String getName() {
    return name;
  }

  public int getSize() {
    return cells.size();
  }

  public int getTreeSize() {
    return rtree.size();
  }

  public int getNodeSize() {
    int i = 0;
    for (Cell c : cells) {
      if (c.hasNodes()) {
        i += c.getNodes().size();
      }
    }
    return i;
  }

  public double[][] getNodeLatLons() {
    double[][] ll = new double[this.getUniqueNodeSize()][2];
    double[] d = new double[2];
    for (Node n : getUniqueNodes()) {
      d[0] = n.getGeoPoint().getLatitude();
      d[1] = n.getGeoPoint().getLongitude();
      ll[ll.length - 1] = d;
    }
    return ll;
  }
  public int[] getNodeIndexes() {
    int[] in  = new int[this.getUniqueNodeSize()];
    for (Node n : getUniqueNodes()) {
      in[in.length - 1] = n.getDataIndex();
    }
    return in;
  }
  public ArrayList<Node> getUniqueNodes() {
    ArrayList<Node> ents = new ArrayList<Node>();
    for (Cell c : cells) {
      if (c.hasNodes()) {
        ents.addAll(c.getNodes());
      }
    }
    Set set = new HashSet(ents);
    ArrayList unique = new ArrayList(set);
    return unique;
  }

  public int getUniqueNodeSize() {
    return getUniqueNodes().size();
  }

  public ArrayList<Edge> getUniqueEdges() {
    ArrayList<Edge> ents = new ArrayList<Edge>();
    for (Cell c : cells) {
      if (c.hasEdges()) {
        ents.addAll(c.getEdges());
      }
    }
    Set set = new HashSet(ents);
    ArrayList unique = new ArrayList(set);
    return unique;
  }

  public int getEdgeSize() {
    int i = 0;
    for (Cell c : cells) {
      if (c.hasEdges()) {
        i += c.getEdges().size();
      }
    }
    return i;
  }

  public ArrayList<Face> getUniqueFaces() {
    ArrayList<Face> ents = new ArrayList<Face>();
    for (Cell c : cells) {
      if (c.hasFaces()) {
        ents.addAll(c.getFaces());
      }
    }
    Set set = new HashSet(ents);
    ArrayList unique = new ArrayList(set);
    return unique;
  }

  public int getFaceSize() {
    int i = 0;
    for (Cell c : cells) {
      if (c.hasFaces()) {
        i += c.getFaces().size();
      }
    }
    return i;
  }

  public LatLonRect getLatLonBoundingBox() {
    LatLonRectangle2D bounds = rtree.getBounds();
    return new LatLonRect((LatLonPoint) new LatLonPointImpl(bounds.getLatMin(), bounds.getLonMin()), (LatLonPoint) new LatLonPointImpl(bounds.getLatMax(), bounds.getLonMax()));
  }

  public Cell getCellFromLatLon(double lat, double lon) {
    LatLonPoint2D p = new LatLonPoint2D.Double(lat, lon);
    return cells.get(rtree.nearest(p));
  }

  public Cell getCellFromLatLon(LatLonPoint2D p) {
    return cells.get(rtree.nearest(p));
  }

  public ArrayList<Cell> getCellsInPolygon(LatLonPolygon2D p) {
    IntArrayList polys = rtree.intersects(p);
    ArrayList<Cell> containedCells = new ArrayList<Cell>(polys.size());
    for (int i : polys.elements()) {
      containedCells.add(this.cells.get(i));
    }
    return containedCells;
  }

  public List<Cell> getCells() {
    return cells;
  }

  public List<Cell> getPolygons() {
    return cells;
  }

  public ConnectivityVariable getConnectivityVariable() {
    return connectivity_variable;
  }

  public Mesh subset(LatLonRect bounds) {
    LatLonRectangle2D r = new LatLonRectangle2D(new LatLonPoint2D.Double(bounds.getUpperLeftPoint().getLatitude(), bounds.getUpperLeftPoint().getLongitude()), new LatLonPoint2D.Double(bounds.getLowerRightPoint().getLatitude(), bounds.getLowerRightPoint().getLongitude()));
    LatLonPolygon2D p = new LatLonPolygon2D.Double(r);
    ArrayList<Cell> containedCells = this.getCellsInPolygon(p);
    return null;
  }

  @Override
  public String toString() {
    // Don't use commas (,) in the string output.
    ArrayList<String> sb = new ArrayList<String>();
    sb.add(this.getName());
    sb.add("Mesh contains: " + this.getSize() + " cells (polygons)");
    sb.add("Mesh contains: " + this.getNodeSize() + " nodes (" + this.getUniqueNodeSize() + " unique) ");
    if (this.getSize() != 0) {
      sb.add("Mesh contains: " + this.getNodeSize() / this.getSize() + " nodes per cell");
    }
    sb.add("Mesh contains: " + this.getEdgeSize() + " edges.");
    sb.add("Mesh contains: " + this.getFaceSize() + " faces.");
    return sb.toString().replace(",","\n");
  }

}
