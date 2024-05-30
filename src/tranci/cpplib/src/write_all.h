// writes all the matrices in a file

write_matrix("sx.op",sx);  // Sx
write_matrix("sy.op",sy);  // Sy
write_matrix("sz.op",sz);  // Sz
write_matrix("lx.op",lx);  // Lx
write_matrix("ly.op",ly);  // Ly
write_matrix("lz.op",lz);  // Lz
write_matrix("ls.op",ls);  // L.S
write_matrix("x2y2.op",x2y2); // x2y2
write_matrix("x2.op",x2); // x2
write_matrix("y2.op",y2); // y2
write_matrix("z2.op",z2); // z2
//write_matrix("xy.op",xy);
//write_matrix("yz.op",yz);
//write_matrix("zx.op",zx);
write_matrix("x4.op",x4); // x4
write_matrix("y4.op",y4); // y4
write_matrix("z4.op",z4); // z4
write_matrix("jx.op",jx); // jx
write_matrix("jy.op",jy);
write_matrix("jz.op",jz);
write_matrix("vc.op",vc);  // coulomb

int iwrite;
stringstream ss;
string sfn;
for (iwrite=0;iwrite<norb;iwrite++)  {
  ss.str(""); // clean stream
  ss << iwrite; // write in stream
  sfn = ss.str(); // convert to string
  write_matrix(sfn+".op",proji[iwrite]); // write in file
  };

