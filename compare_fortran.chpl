use IO;
use parameters;
use domains;
use arrays;

proc load_fortran_grid(fort_file : string, inout arr : [?D] real) {

var tmp_arr : [D3] real;

for k in 0..#nz {
  var f = open("../" + fort_file + ((k+1) : string) + ".dat", iomode.r);
  var r = f.reader(kind=ionative);
  for j in 0..#ny {
    for i in 0..#nx {
      var tmp : real;
      r.readBinary(tmp);
      tmp_arr[k,j,i] = tmp;
    }
  }
  r.close();
}

arr = tmp_arr;

}

proc load_1layer_test(inout arr : [?D] real) {

var tmp_arr : [D] real;

var f = open("../test_grid.dat", iomode.r);
var r = f.reader(kind=ionative);
for j in 0..#ny {
  for i in 0..#nx {
    var tmp : real;
    r.readBinary(tmp);
    tmp_arr[0,j,i] = tmp;
  }
}
r.close();

for k in 1..#nz {
  tmp_arr[k,..,..] = tmp_arr[0,..,..];
}

arr = tmp_arr;

}


proc print_array_3D(arr : [?dom]) {

for k in 0..#nz {
  writeln("layer " + k : string);
  writeln();
  for i in dom[0,..,0] {
    writeln(arr[k,i,..] : real);
    writeln();
  }
}
}

proc print_array_2D(arr : [?dom]) {

for i in dom[..,0] {
  writeln(arr[i,..] : real);
  writeln();
}

}

proc print_array_2D_i(arr : [?dom]) {

for i in dom[..,0] {
  var tmp = arr[i,..];
  tmp = -1i * tmp;
  writeln(tmp : real);
  writeln();
}

}


proc difference3D(fort_file : string, chapel_arr : [?dom]) {

var tmp_arr : [D3] real;

for k in 0..#nz {
  var f = open("../" + fort_file + ((k+1) : string) + ".dat", iomode.r);
  var r = f.reader(kind=ionative);
  for j in 0..#ny {
    for i in 0..#nx {
      var tmp : real;
      r.readBinary(tmp);
      tmp_arr[k,j,i] = tmp;
    }
  }
  r.close();
}

var diff = chapel_arr - tmp_arr;
print_array_3D(diff);
}


proc difference3D_hat(fort_file : string, chapel_arr : [?dom]) {

var tmp_arr : [D3_hat] real;

for k in 0..#nz {
  var f = open("../" + fort_file + ((k+1) : string) + ".dat", iomode.r);
  var r = f.reader(kind=ionative);
  for j in 0..#ny {
    for i in 0..#nx2p {
      var tmp : real;
      r.readBinary(tmp);
      tmp_arr[k,i,j] = tmp;
    }
  }
  r.close();
}

var diff = (chapel_arr : real) - tmp_arr;
print_array_3D(diff);
}


proc difference2D(fort_file : string, chapel_arr : [?dom]) {

var tmp_arr : [D] real;

var f = open("../" + fort_file, iomode.r);
var r = f.reader(kind=ionative);
for j in 0..#ny {
  for i in 0..#nx {
    var tmp : real;
    r.readBinary(tmp);
    tmp_arr[j,i] = tmp;
  }
}

var diff = chapel_arr - tmp_arr;
print_array_2D(diff);
}


proc difference2D_hat(fort_file : string, chapel_arr : [?dom]) {

var tmp_arr : [D_hat] real;

var f = open("../" + fort_file, iomode.r);
var r = f.reader(kind=ionative);
for j in 0..#ny {
  for i in 0..#nx2p {
    var tmp : real;
    r.readBinary(tmp);
    tmp_arr[i,j] = tmp;
  }
}

var diff = chapel_arr - tmp_arr;
print_array_2D(diff);
}


proc difference_nz(fort_file : string, chapel_arr : [?dom]) {

var tmp_arr : [zl2] real;

var f = open("../" + fort_file, iomode.r);
var r = f.reader(kind=ionative);
for j in 0..#nz {
  for i in 0..#nz {
    var tmp : real;
    r.readBinary(tmp);
    tmp_arr[j,i] = tmp;
  }
}

var diff = chapel_arr - tmp_arr;
print_array_2D(diff);
}
