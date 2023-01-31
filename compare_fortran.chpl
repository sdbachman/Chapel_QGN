use IO;
use parameters;
use domains;
use arrays;

proc load_fortran_grid(inout arr : [?D] real) {

var tmp_arr : [D3] real;

var f = open("../../test_grid1.dat", iomode.r);
var r = f.reader(kind=ionative);
for j in 1..ny {
  for i in 1..nx {
    var tmp : real;
    r.readBinary(tmp);
    tmp_arr[1,j,i] = tmp;
  }
}
r.close();

f = open("../../test_grid2.dat", iomode.r);
r = f.reader(kind=ionative);
for j in 1..ny {
  for i in 1..nx {
    var tmp : real;
    r.readBinary(tmp);
    tmp_arr[2,j,i] = tmp;
  }
}
r.close();

f = open("../../test_grid3.dat", iomode.r);
r = f.reader(kind=ionative);
for j in 1..ny {
  for i in 1..nx {
    var tmp : real;
    r.readBinary(tmp);
    tmp_arr[3,j,i] = tmp;
  }
}
r.close();

arr = tmp_arr;

}

proc load_1layer_test(inout arr : [?D] real) {

var tmp_arr : [D] real;

var f = open("../test_grid.dat", iomode.r);
var r = f.reader(kind=ionative);
for j in 1..ny {
  for i in 1..nx {
    var tmp : real;
    r.readBinary(tmp);
    tmp_arr[1,j,i] = tmp;
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

proc print_array_2D_i(arr : [?dom]) {

for i in dom[..,1] {
  var tmp = arr[i,..];
  tmp = -1i * tmp;
  writeln(tmp : real);
  writeln();
}

}


proc difference(fort_file : string, chapel_arr : [?dom]) {

var tmp_arr : [D3] real;

var f = open("../../" + fort_file + "1.dat", iomode.r);
var r = f.reader(kind=ionative);
for j in 1..ny {
  for i in 1..nx {
    var tmp : real;
    r.readBinary(tmp);
    tmp_arr[1,j,i] = tmp;
  }
}
r.close();

f = open("../../" + fort_file + "2.dat", iomode.r);
r = f.reader(kind=ionative);
for j in 1..ny {
  for i in 1..nx {
    var tmp : real;
    r.readBinary(tmp);
    tmp_arr[2,j,i] = tmp;
  }
}
r.close();

f = open("../../" + fort_file + "3.dat", iomode.r);
r = f.reader(kind=ionative);
for j in 1..ny {
  for i in 1..nx {
    var tmp : real;
    r.readBinary(tmp);
    tmp_arr[3,j,i] = tmp;
  }
}
r.close();

var diff = chapel_arr - tmp_arr;
print_array_3D(diff);
}
