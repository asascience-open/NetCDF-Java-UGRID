/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ucar.nc2.dt.ugrid;

import cern.colt.list.IntArrayList;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
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
  private ArrayList<VariableEnhanced> connectivity_variables = new ArrayList<VariableEnhanced>();

  public Mesh(NetcdfDataset ds, VariableEnhanced v) {
    name = v.getName();

    Properties props = new Properties();
    props.setProperty("MaxNodeEntries", "30");
    props.setProperty("MinNodeEntries", "15");
    rtree = new RTree(props);

    processConnectivityVariables(ds, v);
  }

  private void processConnectivityVariables(NetcdfDataset ncd, VariableEnhanced v) {
    ArrayList<VariableEnhanced> foundConns;
    connectivity_variables = null;
    Attribute att = null;
    if (v.findAttributeIgnoreCase("locations") != null) {
      String[] refs = v.findAttributeIgnoreCase("locations").getStringValue().split(" ");
      foundConns = new ArrayList<VariableEnhanced>(refs.length);
      for (String prefix : refs) {
        att = v.findAttributeIgnoreCase(prefix + "_coordinates");
        Attribute conn = v.findAttributeIgnoreCase(prefix + "_connectivity");
        Variable connection_variable;
        if (att != null && conn != null && !att.getStringValue().equals("") && !conn.getStringValue().equals("")) {
          connection_variable = ncd.findVariable(conn.getStringValue());
          foundConns.add((VariableEnhanced)connection_variable);
        } else {
          System.out.println("Could not find coordinate and connectivity information for this Mesh.");
        }
      }
      Set set = new HashSet(foundConns);
      connectivity_variables = new ArrayList(set);
      // A Mesh should only have one connectivity array referenced... right?
      // If not, we may need to generate multiple rtrees.
      constructCells(ncd, att, connectivity_variables.get(0));
    } else {
      System.out.println("No 'locations' attribute on the Mesh.");
    }
  }

  private void constructCells(NetcdfDataset ds, Attribute att, VariableEnhanced connectivity) {
    for (CoordinateSystem coord : ds.getCoordinateSystems()) {
      if (coord.getName().toLowerCase().contains(att.getStringValue().toLowerCase())) {
        // Found the coordinate system
        Attribute cell_type = connectivity.findAttributeIgnoreCase("cell_type");
        if (cell_type != null && !cell_type.getStringValue().equals("")) {
          // for future use of the cell_type attribute

          int start_index;
          if (connectivity.findAttributeIgnoreCase("index_origin") != null) {
            start_index = connectivity.findAttributeIgnoreCase("index_origin").getNumericValue().intValue();
          } else {
            start_index = 0;
            //System.out.println("Connectivity array " + connectivity.getName() + " assumed to start at index 0.");
          }
          if (connectivity.getDimensions().get(0).getLength() > connectivity.getDimensions().get(1).getLength()) {
            //  [N0] [1,2,3]
            //  [N1] [4,5,6]
            //  [N2] [2,5,6]
            //  [N3] [1,4,3]
            processConnectivityArray(connectivity, start_index, coord);
          } else {
            //       N0  N1  N2  N3
            //  --------------------
            //       ̪   ̪   ̪   ̪
            //  [0]  1   4   2   1
            //  [1]  2   5   5   4
            //  [2]  3   6   6   3
            //       ̺   ̺   ̺   ̺
            processConnectivityList(connectivity, start_index, coord);
          }
        } else {
          System.out.println("cell_type was not specified on the "
                  + "connectivity_array variable: " + cell_type.getName());
        }
        break;
      }
    }
  }

  private void processConnectivityList(VariableEnhanced connectivity, int start_index, CoordinateSystem coord) {
    try {
      Node node;
      Cell cell;
      int index;
      ArrayList<Node> nodes = new ArrayList<Node>();
      double[] llats = (double[]) coord.getLatAxis().read().get1DJavaArray(double.class);
      double[] llons = (double[]) coord.getLonAxis().read().get1DJavaArray(double.class);

      ArrayList<Node> unique_nodes = new ArrayList<Node>(llats.length);
      int[][] conn_data = (int[][]) connectivity.read().copyToNDJavaArray();

      for (int j = 0; j < llats.length; j++) {
        node = new Node();
        node.setDataIndex(j);
        node.setGeoPoint(new LatLonPoint2D.Double(llats[j], llons[j]));
        unique_nodes.add(node);
      }
      for (int i = 0; i < connectivity.getDimensions().get(1).getLength(); i++) {
        cell = new Cell();
        nodes.clear();
        for (int k = 0; k < connectivity.getDimensions().get(0).getLength(); k++) {
          index = conn_data[k][i] - start_index;
          nodes.add(unique_nodes.get(index));
        }
        cell.setNodes((ArrayList<Node>) nodes.clone());
        cells.add(cell);
      }
    } catch (IOException e) {
      System.out.println(e);
    }
  }

  private void processConnectivityArray(VariableEnhanced connectivity, int start_index, CoordinateSystem coord) {
    try {
      Node node;
      Cell cell;
      int index;
      ArrayList<Node> nodes = new ArrayList<Node>();
      double[] llats = (double[]) coord.getLatAxis().read().get1DJavaArray(double.class);
      double[] llons = (double[]) coord.getLonAxis().read().get1DJavaArray(double.class);

      ArrayList<Node> unique_nodes = new ArrayList<Node>(llats.length);
      int[][] conn_data = (int[][]) connectivity.read().copyToNDJavaArray();

      for (int j = 0; j < llats.length; j++) {
        node = new Node();
        node.setDataIndex(j);
        node.setGeoPoint(new LatLonPoint2D.Double(llats[j], llons[j]));
        unique_nodes.add(node);
      }

      for (int i = 0; i < connectivity.getDimensions().get(0).getLength(); i++) {
        cell = new Cell();
        nodes.clear();
        for (int k = 0; k < connectivity.getDimensions().get(1).getLength(); k++) {
          index = conn_data[i][k] - start_index;
          nodes.add(unique_nodes.get(index));
        }
        cell.setNodes((ArrayList<Node>) nodes.clone());
        cells.add(cell);
      }

    } catch (IOException e) {
      System.out.println(e);
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

  public List<VariableEnhanced> getConnectivityVariables() {
    return connectivity_variables;
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
