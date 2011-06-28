/*
 * Copyright 1998-2009 University Corporation for Atmospheric Research/Unidata
 *
 * Portions of this software were developed by the Unidata Program at the
 * University Corporation for Atmospheric Research.
 *
 * Access and use of this software shall impose the following obligations
 * and understandings on the user. The user is granted the right, without
 * any fee or cost, to use, copy, modify, alter, enhance and distribute
 * this software, and any derivative works thereof, and its supporting
 * documentation for any purpose whatsoever, provided that this entire
 * notice appears in all copies of the software, derivative works and
 * supporting documentation.  Further, UCAR requests that the user credit
 * UCAR/Unidata in any publications that result from the use of this
 * software or in any product that includes this software. The names UCAR
 * and/or Unidata, however, may not be used in any advertising or publicity
 * to endorse or promote any products or commercial entity unless specific
 * written permission is obtained from UCAR/Unidata. The user also
 * understands that UCAR/Unidata is not obligated to provide the user with
 * any support, consulting, training or assistance of any kind with regard
 * to the use, operation and performance of this software nor to provide
 * the user with any updates, revisions, new versions or "bug fixes."
 *
 * THIS SOFTWARE IS PROVIDED BY UCAR/UNIDATA "AS IS" AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL UCAR/UNIDATA BE LIABLE FOR ANY SPECIAL,
 * INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING
 * FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
 * NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION
 * WITH THE ACCESS, USE OR PERFORMANCE OF THIS SOFTWARE.
 */
package ucar.nc2.dt.ugrid;

import ucar.nc2.dataset.*;
import ucar.nc2.Attribute;
import ucar.nc2.VariableSimpleIF;
import ucar.nc2.NetcdfFile;
import ucar.nc2.Variable;
import ucar.nc2.dt.UGridDatatype;
import ucar.nc2.util.cache.FileCache;
import ucar.nc2.constants.FeatureType;
import ucar.nc2.units.DateRange;
import ucar.unidata.geoloc.LatLonRect;

import java.util.*;
import java.io.IOException;
import ucar.nc2.constants.CF;

/**
 * Make a NetcdfDataset into a collection of GeoGrids with Georeferencing coordinate systems.
 * <p/>
 * <p/>
 * A variable will be made into a GeoGrid if it has a Georeferencing coordinate system,
 * using GridCoordSys.isGridCoordSys(), and it has no extra dimensions, ie
 * GridCoordSys.isComplete( var) is true.
 * If it has multiple Georeferencing coordinate systems, any one that is a product set will be given preference.
 * <p/>
 * Example:
 * <p/>
 * <pre>
 * GridDataset gridDs = GridDataset.open (uriString);
 * List grids = gridDs.getGrids();
 * for (int i=0; i&lt;grids.size(); i++) {
 *   GeoGrid grid = (Geogrid) grids.get(i);
 * }
 * </pre>
 *
 * @author caron
 */

public class UGridDataset implements ucar.nc2.dt.UGridDataset, ucar.nc2.ft.FeatureDataset {
  private NetcdfDataset ds;
  private ArrayList<MeshVariable> meshVariables = new ArrayList<MeshVariable>();
  private Map<String, Meshset> meshsetHash = new HashMap<String, Meshset>();
  
  /**
   * Open a netcdf dataset, using NetcdfDataset.defaultEnhanceMode plus CoordSystems
   * and turn into a UGridDataset.
   *
   * @param location netcdf dataset to open, using NetcdfDataset.acquireDataset().
   * @return GridDataset
   * @throws java.io.IOException on read error
   * @see ucar.nc2.dataset.NetcdfDataset#acquireDataset
   */
  static public UGridDataset open(String location) throws java.io.IOException {
    return open(location, NetcdfDataset.getDefaultEnhanceMode());
  }

  /**
   * Open a netcdf dataset, using NetcdfDataset.defaultEnhanceMode plus CoordSystems
   * and turn into a UGridDataset.
   *
   * @param location netcdf dataset to open, using NetcdfDataset.acquireDataset().
   * @param enhanceMode open netcdf dataset with this enhanceMode
   * @return GridDataset
   * @throws java.io.IOException on read error
   * @see ucar.nc2.dataset.NetcdfDataset#acquireDataset
   */
  static public UGridDataset open(String location, Set<NetcdfDataset.Enhance> enhanceMode) throws java.io.IOException {
    NetcdfDataset ds = ucar.nc2.dataset.NetcdfDataset.acquireDataset(null, location, enhanceMode, -1, null, null);
    return new UGridDataset(ds);
  }

  /**
   * Create a UGridDataset from a NetcdfDataset.
   *
   * @param ds underlying NetcdfDataset, will do Enhance.CoordSystems if not already done.
   * @throws java.io.IOException on read error
   */
  public UGridDataset(NetcdfDataset ds) throws IOException {
    this(ds, null);
  }

  /**
   * Create a UGridDataset from a NetcdfDataset.
   *
   * @param ds underlying NetcdfDataset, will do Enhance.CoordSystems if not already done.
   * @param parseInfo put parse info here, may be null
   * @throws java.io.IOException on read error
   */
  public UGridDataset(NetcdfDataset ds, Formatter parseInfo) throws IOException {
    this.ds = ds;
    ds.enhance(NetcdfDataset.getDefaultEnhanceMode());
    // look for Meshes
    if (parseInfo != null) parseInfo.format("UGridDataset looking for MeshVariables\n");
    List<Variable> vars = ds.getVariables();
    for (Variable var : vars) {
      // See how many "Mesh" are defined in the dataset
      // "Mesh" is defined by the standard_name "topology_description"
      if ((var.findAttributeIgnoreCase(CF.STANDARD_NAME)) != null && (var.findAttributeIgnoreCase(CF.STANDARD_NAME).getStringValue().equals("topology_description"))) {
        VariableEnhanced varDS = (VariableEnhanced) var;
        constructMeshVariable(ds, varDS, parseInfo);
      }
    }
  }

  private void constructMeshVariable(NetcdfDataset ds, VariableEnhanced v, Formatter parseInfo) {
    // Add Mesh to the "meshes" Hash
    if (v instanceof StructureDS) {
      StructureDS s = (StructureDS) v;
      List<Variable> members = s.getVariables();
      for (Variable nested : members) {
        constructMeshVariable(ds, (VariableEnhanced) nested, parseInfo);
      }
    } else {
      Mesh m = new Mesh(ds,v);
      if (m != null) {
        addMesh((VariableDS) v, m, parseInfo);
      }
    }
  }

  private void addMesh(VariableDS varDS, Mesh m, Formatter parseInfo) {
    Meshset meshset;
    if (null == (meshset = meshsetHash.get(m.getName()))) {
      meshset = new Meshset(m, varDS);
      meshsetHash.put(m.getName(), meshset);
      if (parseInfo != null) parseInfo.format(" -make new Mesh= %s\n", m.getName());
      //m.makeVerticalTransform(this, parseInfo);
      setVariables(meshset);
    }
  }

  private void setVariables(Meshset meshset) {
    for (Variable v : ds.getVariables()) {
      if (v.findAttributeIgnoreCase("mesh") != null) {
        if (v.findAttributeIgnoreCase("mesh").getStringValue().toLowerCase().contains(meshset.getMesh().getName().toLowerCase())) {
          MeshVariable mv = new MeshVariable(this, (VariableDS)v, meshset);
          if (!meshVariables.contains(mv)) {
            meshVariables.add(mv);
            meshset.add(mv);
          } else {
            mv = null;
          }
        }
      }
    }
  }
  
  public UGridDatatype getMeshVariableByName(String name) {
    UGridDatatype z = null;
    for (UGridDatatype m : meshVariables) {
      if (m.getName().equals(name)) {
        z = m;
      }
    }
    return z;
  }

  public List<UGridDatatype> getMeshVariables() {
    return new ArrayList<UGridDatatype>(meshVariables);
  }

  public void calcBounds() throws java.io.IOException {
    // not needed
  }

  public List<Attribute> getGlobalAttributes() {
    return ds.getGlobalAttributes();
  }

  public String getTitle() {
    String title = ds.findAttValueIgnoreCase(null, "title", null);
    return (title == null) ? getName() : title;
  }

  public String getDescription() {
    String desc = ds.findAttValueIgnoreCase(null, "description", null);
    if (desc == null)
      desc = ds.findAttValueIgnoreCase(null, "history", null);
    return (desc == null) ? getName() : desc;
  }

  public String getName() {
    return ds.getLocation();
  }

  public String getLocation() {
    return ds.getLocation();
  }

  public String getLocationURI() {
    return ds.getLocation();
  }

  private DateRange dateRangeMax = null;
  private LatLonRect llbbMax = null;

  private void makeRanges() {

    for (ucar.nc2.dt.UGridDataset.Meshset ms : getMeshsets()) {
      Mesh m = ms.getMesh();
      LatLonRect llbb = m.getLatLonBoundingBox();
      if (llbbMax == null) {
        llbbMax = llbb;
      } else {
        llbbMax.extend(llbb);
      }
      /*
      DateRange dateRange = m.getDateRange();
      if (dateRange != null) {
        if (dateRangeMax == null)
          dateRangeMax = dateRange;
        else
          dateRangeMax.extend(dateRange);
      }
      */
    }
  }

  public Date getStartDate() {
    if (dateRangeMax == null) makeRanges();
    return (dateRangeMax == null) ? null : dateRangeMax.getStart().getDate();
  }

  public Date getEndDate() {
    if (dateRangeMax == null) makeRanges();
    return (dateRangeMax == null) ? null : dateRangeMax.getEnd().getDate();
  }

  public LatLonRect getBoundingBox() {
    if (llbbMax == null) makeRanges();
    return llbbMax;
  }

  public Attribute findGlobalAttributeIgnoreCase(String name) {
    return ds.findGlobalAttributeIgnoreCase(name);
  }

  public List<VariableSimpleIF> getDataVariables() {
    List<VariableSimpleIF> result = new ArrayList<VariableSimpleIF>(meshVariables.size());
    for (UGridDatatype mv : getMeshVariables()) {
      if (mv.getVariable() != null)
        result.add(mv.getVariable());
    }
    return result;
  }

  public VariableSimpleIF getDataVariable(String shortName) {
    return ds.findTopVariable(shortName);
  }

  public NetcdfFile getNetcdfFile() {
    return ds;
  }

  public FeatureType getFeatureType() {
    return FeatureType.UGRID;
  }

  public DateRange getDateRange() {
    if (dateRangeMax == null) makeRanges();
    return dateRangeMax;
  }

  public String getImplementationName() {
    return ds.getConventionUsed();
  }


  public synchronized void close() throws java.io.IOException {
    if (fileCache != null) {
      fileCache.release(this);
    } else {
      try {
        if (ds != null) ds.close();
      } finally {
        ds = null;
      }
    }
  }

  public boolean sync() throws IOException {
    return (ds != null) ? ds.sync() : false;
  }

  protected FileCache fileCache;
  public void setFileCache(FileCache fileCache) {
    this.fileCache = fileCache;
  }

  // TODO: Show info about Grid and Coord systems
  private void getInfo(Formatter buf) {
  }

  public String getDetailInfo() {
    Formatter buff = new Formatter();
    getDetailInfo(buff);
    return buff.toString();
  }

  public NetcdfDataset getNetcdfDataset() {
    return ds;
  }

  public void getDetailInfo(Formatter buff) {
    getInfo(buff);
    buff.format("\n\n----------------------------------------------------\n");
    NetcdfDatasetInfo info = null;
    try {
      info = new NetcdfDatasetInfo( ds.getLocation());
      buff.format("%s", info.getParseInfo());
    } catch (IOException e) {
      buff.format("NetcdfDatasetInfo failed");
    } finally {
      if (info != null) try { info.close(); } catch (IOException ee) {} // do nothing
    }
    buff.format("\n\n----------------------------------------------------\n");
    buff.format("%s", ds.toString());
    buff.format("\n\n----------------------------------------------------\n");
  }

  public List<ucar.nc2.dt.UGridDataset.Meshset> getMeshsets() {
    return new ArrayList<ucar.nc2.dt.UGridDataset.Meshset>(meshsetHash.values());
  }

  /**
   * This is a set of MeshVariables with the same Mesh and Topology
   */
  public class Meshset implements ucar.nc2.dt.UGridDataset.Meshset {

    private Mesh mesh;
    private VariableDS description_variable;
    private List<UGridDatatype> meshVariables = new ArrayList<UGridDatatype>();

    public Meshset(Mesh m, VariableDS conn) {
      this.mesh = m;
      this.description_variable = conn;
    }

    private void add(MeshVariable mv) {
      meshVariables.add(mv);
    }

    /**
     * Get list of MeshVariable objects
     */
    public List<UGridDatatype> getMeshVariables() {
      return meshVariables;
    }

    public UGridDatatype getMeshVariableByName(String name) {
      UGridDatatype z = null;
      for (UGridDatatype m : meshVariables) {
        if (m.getName().equals(name)) {
          z = m;
        }
      }
      return z;
    }

    public Mesh getMesh() {
      return mesh;
    }

    public VariableDS getDescriptionVariable() {
      return description_variable;
    }

  }
}