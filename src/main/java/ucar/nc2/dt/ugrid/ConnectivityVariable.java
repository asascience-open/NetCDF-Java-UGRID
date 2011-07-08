/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ucar.nc2.dt.ugrid;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import ucar.ma2.Array;
import ucar.ma2.DataType;
import ucar.nc2.Attribute;
import ucar.nc2.Dimension;
import ucar.nc2.Variable;
import ucar.nc2.constants.AxisType;
import ucar.nc2.dataset.CoordinateSystem;
import ucar.nc2.dataset.NetcdfDataset;
import ucar.nc2.dataset.VariableDS;
import ucar.nc2.dt.ugrid.geom.LatLonPoint2D;

/**
 *
 * @author Kyle
 */
public class ConnectivityVariable {
  
  private enum Type {WIDE,TALL}
  private final static String CELL_TYPE = "cell_type";
  private final static String INDEX_ORIGIN = "index_origin";
  
  private Type type;
  private String cellType = "none";
  private int startIndex = 0;
  private Variable variable;

  public ConnectivityVariable(Variable var) {
    variable = var;
    parse();
  }
  
  private void parse() {
    /*
     * Cell Type
     * For future use of the cell_type attribute
     * "tri_ccw" seems to map to triangles!
     */
    Attribute cell_type = variable.findAttributeIgnoreCase(CELL_TYPE);
    if (cell_type != null && !cell_type.getStringValue().isEmpty()) {
      cellType = cell_type.getStringValue();
    }
    
    /*
     * Start Index
     */
    if (variable.findAttributeIgnoreCase(INDEX_ORIGIN) != null) {
      startIndex = variable.findAttributeIgnoreCase(INDEX_ORIGIN).getNumericValue().intValue();
    }
    
    /*
     * Connectivity Array
     * There are two types of ways to specify a connectivity array!!
     */
    if (variable.getDimensions().get(0).getLength() > variable.getDimensions().get(1).getLength()) {
      //  [N0] [1,2,3]
      //  [N1] [4,5,6]
      //  [N2] [2,5,6]
      //  [N3] [1,4,3]
      type = Type.TALL;
    } else {
      //       N0  N1  N2  N3
      //  --------------------
      //       ̪   ̪   ̪   ̪
      //  [0]  1   4   2   1
      //  [1]  2   5   5   4
      //  [2]  3   6   6   3
      //       ̺   ̺   ̺   ̺
      type = Type.WIDE;
    }
  }
  
  public List<Cell> createCells(List<String> locations, List<CoordinateSystem> coords) {
    try {
      CoordinateSystem node_cs = null;
      CoordinateSystem face_cs = null;

      for (int i = 0 ; i < locations.size() ; i++) {
        if (locations.get(i).equalsIgnoreCase("node")) {
          node_cs = coords.get(i);
        } else if (locations.get(i).equalsIgnoreCase("face")) {
          face_cs = coords.get(i);
        } 
//        else if (locations.get(i).equalsIgnoreCase("edge")) {
//          edge_cs = coords.get(i);
//        }
      }

      List<Cell> cells = new ArrayList<Cell>();
      Node node;
      Face face;
      Edge edge;
      Cell cell;
      int index;
      int[][] conn_data = (int[][]) variable.read().copyToNDJavaArray();
      
      double[] face_lats = null;
      double[] face_lons = null;
      ArrayList<Face> unique_faces = null;
      if (face_cs != null) {
        face_lats = (double[]) face_cs.getLatAxis().read().get1DJavaArray(double.class);
        face_lons = (double[]) face_cs.getLonAxis().read().get1DJavaArray(double.class);
        unique_faces = new ArrayList<Face>(face_lats.length);
        for (int g = 0; g < face_lats.length; g++) {
          face = new Face();
          face.setDataIndex(g);
          face.setGeoPoint(new LatLonPoint2D.Double(face_lats[g], face_lons[g]));
          unique_faces.add(face);
        }
      }
      
      double[] node_lats = (double[]) node_cs.getLatAxis().read().get1DJavaArray(double.class);
      double[] node_lons = (double[]) node_cs.getLonAxis().read().get1DJavaArray(double.class);
      ArrayList<Node> unique_nodes = new ArrayList<Node>(node_lats.length);
      for (int j = 0; j < node_lats.length; j++) {
        node = new Node();
        node.setDataIndex(j);
        node.setGeoPoint(new LatLonPoint2D.Double(node_lats[j], node_lons[j]));
        unique_nodes.add(node);
      }
      
      ArrayList<Node> nodes = new ArrayList<Node>();
      ArrayList<Face> faces = new ArrayList<Face>();
      if (this.isTall()) {
        for (int i = 0; i < conn_data.length; i++) {
          cell = new Cell();
          nodes.clear();
          faces.clear();
          for (int k = 0; k < conn_data[i].length; k++) {
            index = conn_data[i][k] - startIndex;
            nodes.add(unique_nodes.get(index));
          }
          if (!cell.hasFaces() && unique_faces != null) {
            faces.add(unique_faces.get(i));
            cell.setFaces(faces);
          }
          cell.setNodes((ArrayList<Node>) nodes.clone());
          cell.setConnectivityIndex(i);
          cells.add(cell);
        }
      } else {
        for (int i = 0; i < conn_data[0].length; i++) {
          cell = new Cell();
          nodes.clear();
          faces.clear();
          for (int k = 0; k < conn_data.length; k++) {
            index = conn_data[k][i] - startIndex;
            nodes.add(unique_nodes.get(index));
          }
          if (!cell.hasFaces() && unique_faces != null) {
            faces.add(unique_faces.get(i));
            cell.setFaces(faces);
          }
          cell.setNodes((ArrayList<Node>) nodes.clone());
          cell.setConnectivityIndex(i);
          cells.add(cell);
        }
      }
      return cells;
    } catch (IOException e) {
      return null;
    }
  }
  
  public void subsetToDataset(UGridDataset ugd, NetcdfDataset ncd, List<Cell> containedCells) {
  
    /*
     * The 'NFaces' dimension... or the number of sides of the polygons in the grid.
     * We are assuming this is 3 in other parts of the code, but not here.
     */
    String cell_number_dimension_name;
    if (this.isTall()) {
      ncd.addDimension(null, this.variable.getDimension(1));
      cell_number_dimension_name = variable.getDimension(0).getName();
    } else {
      ncd.addDimension(null, this.variable.getDimension(0));
      cell_number_dimension_name = variable.getDimension(1).getName();
    }
    ncd.finish();
    
    // Create cell_dim of correct size if it does not exist
    Dimension cell_dim = ncd.findDimension(cell_number_dimension_name);
    if (cell_dim == null) {
      cell_dim = ncd.addDimension(null, new Dimension(cell_number_dimension_name, containedCells.size()));
    }
    ncd.finish();
    
    Variable newConn = new VariableDS(ncd, null, null, variable.getShortName(), DataType.INT, variable.getDimensionsString(), null, null);
    for (Attribute a : (List<Attribute>) variable.getAttributes()) {
      newConn.addAttribute(a);
    }
    ncd.finish();
    
    int[][] raw_data = new int[newConn.getDimension(0).getLength()][newConn.getDimension(1).getLength()];
    for (int i = 0 ; i < containedCells.size() ; i++) {
      for (int j = 0 ; j < containedCells.get(i).getNodes().size() ; j++) {
        if (this.isTall()) {
          raw_data[i][j] = containedCells.get(i).getNodes().get(j).getDataIndex();
        } else {
          raw_data[j][i] = containedCells.get(i).getNodes().get(j).getDataIndex();
        }
      }
    }
    Array conn_data = Array.factory(raw_data);
    newConn.setCachedData(conn_data);
    
    ncd.addVariable(null, newConn);
    ncd.finish();
  }
  
  public VariableDS subsetToVariable(List<Cell> containedCells) {
    Variable newV = new Variable(variable);
    int[][] conn = new int[newV.getShape(0)][newV.getShape(1)];
    int count = 0;
    for (Cell c : containedCells) {
      for (int i = 0 ; i < c.getNodes().size() ; i++) {
        if (this.isTall()) {
          conn[count][i] = c.getNodes().get(i).getDataIndex();
        } else {
          conn[i][count] = c.getNodes().get(i).getDataIndex();
        }
      }
      count++;
    }

    newV.setCachedData(Array.factory(conn));

    return new VariableDS(null, newV, false);
  }
  
//  public HashMap generateIndexMap(List<Cell> containedCells) {
//    HashMap<Integer,Integer> news = new HashMap<Integer,Integer>();
//    for (Cell c : containedCells) {
//      for (int i = 0 ; i < c.getNodes().size() ; i++) {
//        if (this.isTall()) {
//          conn[count][i] = c.getNodes().get(i).getDataIndex();
//        } else {
//          conn[i][count] = c.getNodes().get(i).getDataIndex();
//        }
//      }
//      count++;
//    }
//    return news;
//  }
  
  private Variable getVariable() {
    return variable;
  }

  public void setVariable(Variable variable) {
    this.variable = variable;
    parse();
  }
  
  public boolean isTall() {
    return type.equals(Type.TALL);
  }
  
}
