use IO;
use parameters;
use domains;
use arrays;

proc load_fortran_grid(out arr : [] real) {

var f = open("../../test_grid1.dat", iomode.r);
var r = f.reader(kind=ionative);
for i in 1..nx {
  for j in 1..ny {
      var tmp : real;
      r.readBinary(tmp);
    arr[1,i,j] = tmp;
  }
}
r.close();

f = open("../../test_grid2.dat", iomode.r);
r = f.reader(kind=ionative);
for i in 1..nx {
  for j in 1..ny {
      var tmp : real;
      r.readBinary(tmp);
    arr[1,i,j] = tmp;
  }
}
r.close();

f = open("../../test_grid3.dat", iomode.r);
r = f.reader(kind=ionative);
for i in 1..nx {
  for j in 1..ny {
      var tmp : real;
      r.readBinary(tmp);
    arr[1,i,j] = tmp;
  }
}
r.close();

}

proc load_1layer_test(inout arr : [?D] real) {

var tmp_arr : [D] real;

var f = open("../test_grid.dat", iomode.r);
var r = f.reader(kind=ionative);
for i in 1..nx {
  for j in 1..ny {
    var tmp : real;
    r.readBinary(tmp);
    tmp_arr[1,i,j] = tmp;
  }
}
r.close();

for k in 2..nz {
  tmp_arr[k,..,..] = tmp_arr[1,..,..];
}

arr = tmp_arr;

}


proc print_array_3D(arr : [?dom]) {

for k in 1..nz {
  writeln("layer " + k : string);
  writeln();
  for i in dom[1,..,1] {
    writeln(arr[k,i,..] : real);
    writeln();
  }
}
}

proc print_array_2D(arr : [?dom]) {

for i in dom[..,1] {
  writeln(arr[i,..] : real);
  writeln();
}

}
