use parameters;
use domains;
use arrays;
use diag_table;
use QGN_Module;
use FFT_utils;

use FileSystem;
use BlockDist;
use NetCDF.C_NetCDF;


proc Diagnostics(i : int) {

  if (Q_DIAG || SAVE_CHECKPOINTS) {
    if ( ((i % Q_DIAG_FREQ)==0) || ((i % CHECKPOINT_FREQ)==0) ) {
      var q_tmp : [D3] real;
      execute_backward_FFTs(q_hat, q_tmp);
      normalize(q_tmp);
      WriteOutput(q_tmp, "q", "seconds-1", i);
    }
  }

  if (PSI_DIAG) || (U_DIAG) || (V_DIAG) {
    if ((i % PSI_DIAG_FREQ)==0) || ((i % U_DIAG_FREQ)==0) || ((i % V_DIAG_FREQ)==0) {
      GetPsi(q_hat);
    }
  }

  if (PSI_DIAG) {
    if ((i % PSI_DIAG_FREQ)==0) {
      var psi_tmp : [D3] real;
      execute_backward_FFTs(psi_hat, psi_tmp);
      normalize(psi_tmp);
      WriteOutput(psi_tmp, "psi", "meters2 seconds-1", i);
    }
  }

  if (U_DIAG) {
    if ((i % U_DIAG_FREQ)==0) {
      forall (i,j,k) in D3_hat {
        u_hat[i,j,k] = -1i*ky[j,k]*psi_hat[i,j,k];
      }
      var u_tmp : [D3] real;
      execute_backward_FFTs(u_hat, u_tmp);
      normalize(u_tmp);
      WriteOutput(u_tmp, "u", "meters second-1", i);
    }
  }

  if (V_DIAG) {
    if ((i % V_DIAG_FREQ)==0) {
      forall (i,j,k) in D3_hat {
        v_hat[i,j,k] = 1i*kx[j,k]*psi_hat[i,j,k];
      }
      var v_tmp : [D3] real;
      execute_backward_FFTs(v_hat, v_tmp);
      normalize(v_tmp);
      WriteOutput(v_tmp, "v", "meters second-1", i);
    }
  }

}


proc WriteOutput(ref arr_in: [?D] real, varName : string, units : string, i : int) {

  /* IDs for the netCDF file, dimensions, and variables. */
    var ncid, x_dimid, y_dimid, z_dimid, time_dimid : c_int;
    var x_varid, y_varid, z_varid, time_varid : c_int;
    var varid : c_int;

    var ndims : int = 4;
    var dimids: [0..#ndims] c_int;

    var shape = arr_in.shape;

    var att_text : string;

    var current_time = t;
    var zo = z;
    var yo = y[..,0];
    var xo = x[0,..];

    var timeName = "time";
    var zName = "z";
    var yName = "y";
    var xName = "x";

  /* The timestamp for the filename */
    var currentIter = i : string;
    const maxLen = 10;
    const zero_len = maxLen - currentIter.size;
    const paddedStr = (zero_len * "0") + currentIter;

  /* Create the file. */
    extern proc nc_create(path : c_string, cmode : c_int, ncidp : c_ptr(c_int)) : c_int;
    nc_create( (varName + "." + paddedStr + ".nc").c_str(), NC_CLOBBER, c_ptrTo(ncid));

  /* Define the dimensions. The record dimension is defined to have
     unlimited length - it can grow as needed. In this example it is
     the time dimension.*/
    extern proc nc_def_dim(ncid : c_int, name : c_string, len : c_size_t, idp : c_ptr(c_int)) : c_int;
    nc_def_dim(ncid, timeName.c_str(), 1 : c_size_t, time_dimid);
    nc_def_dim(ncid, zName.c_str(), shape[0] : c_size_t, z_dimid);
    nc_def_dim(ncid, yName.c_str(), shape[1] : c_size_t, y_dimid);
    nc_def_dim(ncid, xName.c_str(), shape[2] : c_size_t, x_dimid);

  /* Define the coordinate variables. */
    extern proc nc_def_var(ncid : c_int, name : c_string, xtype : nc_type, ndims : c_int, dimidsp : c_ptr(c_int), varidp : c_ptr(c_int)) : c_int;
      nc_def_var(ncid, timeName.c_str(), NC_DOUBLE, 1 : c_int, c_ptrTo(time_dimid), c_ptrTo(time_varid));
      nc_def_var(ncid, zName.c_str(), NC_DOUBLE, 1 : c_int, c_ptrTo(z_dimid), c_ptrTo(z_varid));
      nc_def_var(ncid, yName.c_str(), NC_DOUBLE, 1 : c_int, c_ptrTo(y_dimid), c_ptrTo(y_varid));
      nc_def_var(ncid, xName.c_str(), NC_DOUBLE, 1 : c_int, c_ptrTo(x_dimid), c_ptrTo(x_varid));

  /* Assign units attributes to coordinate variables. */
    extern proc nc_put_att_text(ncid : c_int, varid : c_int, name : c_string, len : c_size_t, op : c_string) : c_int;
    att_text = "timesteps";
    nc_put_att_text(ncid, time_varid, "units".c_str(), att_text.numBytes : c_size_t, att_text.c_str());
    att_text = "meters";
    nc_put_att_text(ncid, z_varid, "units".c_str(), att_text.numBytes : c_size_t, att_text.c_str());
    att_text = "meters";
    nc_put_att_text(ncid, y_varid, "units".c_str(), att_text.numBytes : c_size_t, att_text.c_str());
    att_text = "meters";
    nc_put_att_text(ncid, x_varid, "units".c_str(), att_text.numBytes : c_size_t, att_text.c_str());

  /* The dimids array is used to pass the dimids of the dimensions of
     the netCDF variables. In C, the unlimited dimension must come first on the list of dimids. */
    dimids[0] = time_dimid;
    dimids[1] = z_dimid;
    dimids[2] = y_dimid;
    dimids[3] = x_dimid;

  /* Define the netCDF variable. */
    nc_def_var(ncid, varName.c_str(), NC_DOUBLE, ndims : c_int, c_ptrTo(dimids[0]), c_ptrTo(varid));

  /* Assign units attributes to the netCDF variables. */
    nc_put_att_text(ncid, varid, "units".c_str(), units.numBytes : c_size_t, units.c_str());

  /* End define mode. */
    nc_enddef(ncid);

  /* Write the coordinate variable data. */
    extern proc nc_put_var_double(ncid : c_int, varid : c_int, op : c_ptr(c_double)) : c_int;
    nc_put_var_double(ncid, time_varid, c_ptrTo(current_time));
    nc_put_var_double(ncid, z_varid, c_ptrTo(zo[0]));
    nc_put_var_double(ncid, y_varid, c_ptrTo(yo[0]));
    nc_put_var_double(ncid, x_varid, c_ptrTo(xo[0]));

  /* Create arrays of c_size_t for compatibility with NetCDF-C functions. */
  /* Determine where to start reading file, and how many elements to read */
    // Start specifies a hyperslab.  It expects an array of dimension sizes
      var start = tuplify(D.localSubdomain().first);
    // Count specifies a hyperslab.  It expects an array of dimension sizes
      var count = tuplify(D.localSubdomain().shape);

  /* Adding an extra first element to account for the "time" dimension. */
    var start_c : [0..start.size] c_size_t;
    var count_c : [0..count.size] c_size_t;
    start_c[0] = 0 : c_size_t;
    count_c[0] = 1 : c_size_t;
    for i in 1..start.size {
      start_c[i] = start[i-1] : c_size_t;
      count_c[i] = count[i-1] : c_size_t;
    }

    extern proc nc_put_vara_double(ncid : c_int, varid : c_int, startp : c_ptr(c_size_t), countp : c_ptr(c_size_t), op : c_ptr(c_double)) : c_int;
    nc_put_vara_double(ncid, varid, c_ptrTo(start_c), c_ptrTo(count_c), c_ptrTo(arr_in[start]));

    nc_close(ncid);

}

inline proc tuplify(x) {
  if isTuple(x) then return x; else return (x,);
}

proc read_initial_state(ref q_in : [] real) {

  var input_files = glob("q.0*");

  var filename = input_files.last;
  writeln("Reading restart from file ", filename, ".");
  var varName = "q";

  var ncid : c_int;
  var varid : c_int;
  var ndims : c_int = 4;
  var dimid: c_int;

  // Open the file
  // (1)  int nc_open(const char* path, int mode,     int* ncidp)
  extern proc nc_open(path : c_string, mode : c_int, ncidp : c_ptr(c_int)) : c_int;
  nc_open(filename.c_str(), NC_NOWRITE, c_ptrTo(ncid));

  // Get the variable ID
  //
  //      int nc_inq_varid(int ncid,    const char* name,      int* varidp)
  extern proc nc_inq_varid(ncid: c_int, varName: c_string, varid: c_ptr(c_int));
  nc_inq_varid(ncid, varName.c_str(), c_ptrTo(varid));

  var dimids : [0..#ndims] c_int;

  // Get the IDs of each dimension
  //
  //      int nc_inq_vardimid(int ncid,     int varid,     int* dimidsp)
  extern proc nc_inq_vardimid(ncid : c_int, varid : c_int, dimidsp : c_ptr(c_int)) : c_int;

  nc_inq_vardimid(ncid, varid, c_ptrTo(dimids));

  var dimlens : [0..#(ndims-1)] c_size_t;

  // Get the size of each dimension
  //
  //      int nc_inq_dimlen(int ncid,     int dimid,     size_t* lenp)
  extern proc nc_inq_dimlen(ncid : c_int, dimid : c_int, lenp : c_ptr(c_size_t)) : c_int;
  for i in 0..#(ndims-1) do {
    nc_inq_dimlen(ncid, dimids[i+1], c_ptrTo(dimlens[i]));
  }

  param p = 3;
  var dom_in = CreateDomain(p, dimlens);
  q_in = DistributedRead(filename, varid, dom_in);

/* Get the time and timestep from the restart file */
  nc_inq_varid(ncid, "time".c_str(), c_ptrTo(varid));

  var start_c : c_size_t = 0;
  var count_c : c_size_t = 1;
  extern proc nc_get_var1_double(ncid : c_int, varid : c_int, indexp : c_ptr(c_size_t), ip : c_ptr(real(64)));
  nc_get_var1_double(ncid, varid, c_ptrTo(start_c), c_ptrTo(t));
  writeln("Restart time is ", t, ".");

  var fn_split = filename.split(".");
  Nt_start = (fn_split[1] : int);
  writeln("Restart timestep is ", Nt_start, ".");

}


proc DistributedRead(const filename, varid, dom_in) {

  const D = dom_in dmapped Block (dom_in);
  var dist_array : [D] real(64);

  coforall loc in Locales do on loc {
    writeln("Local subdomain on Locale ", here.id, ": \n", D.localSubdomain());

    /* Some external procedure declarations */
      extern proc nc_get_vara_double(ncid : c_int, varid : c_int, startp : c_ptr(c_size_t), countp : c_ptr(c_size_t), ip : c_ptr(real(64))) : c_int;

    /* Determine where to start reading file, and how many elements to read */
      // Start specifies a hyperslab.  It expects an array of dimension sizes
      var start = tuplify(D.localSubdomain().first);
      // Count specifies a hyperslab.  It expects an array of dimension sizes
      var count = tuplify(D.localSubdomain().shape);

    /* Create arrays of c_size_t for compatibility with NetCDF-C functions. */
      var start_c : [0..start.size] c_size_t;
      var count_c : [0..count.size] c_size_t;
      start_c[0] = 0 : c_size_t;
      count_c[0] = 1 : c_size_t;
      for i in 1..start.size {
        start_c[i] = start[i-1] : c_size_t;
        count_c[i] = count[i-1] : c_size_t;
      }

      var ncid : c_int;
      nc_open(filename.c_str(), NC_NOWRITE, ncid);

      writeln("URL on Locale ", here.id, ": ", filename);

      nc_get_vara_double(ncid, varid, c_ptrTo(start_c), c_ptrTo(count_c), c_ptrTo(dist_array[start]));

      nc_close(ncid);
  }

  return dist_array;

}

proc CreateDomain(param numDims, indicesArr) {
  var indices: numDims*range;
  for param i in 0..<numDims do
    indices[i] = 0..#indicesArr[i];
  return {(...indices)};
}
