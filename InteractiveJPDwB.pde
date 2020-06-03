// Interactive demo of the persistence pipeline in 2D,
// building on javaplexDemo.pde from https://github.com/appliedtopology/javaplex.


import edu.stanford.math.plex4.api.*;
import edu.stanford.math.plex4.examples.*;
import edu.stanford.math.plex4.streams.impl.VietorisRipsStream;
import edu.stanford.math.plex4.homology.chain_basis.Simplex;
import edu.stanford.math.plex4.homology.filtration.FiltrationConverter;
import edu.stanford.math.plex4.homology.interfaces.AbstractPersistenceAlgorithm;
import edu.stanford.math.plex4.homology.barcodes.*;
import java.util.Map.Entry;
import java.util.List;

double[][] pts;
double[][] pts1;

float offsetX,offsetY,sizeX,sizeY;
int dragX,dragY,oldmouseX,oldmouseY;
double eps = 0.01;
double f = eps;
double maxeps = 0.3;
Table table;
VietorisRipsStream<double[]> vrs;
FiltrationConverter fc;
AbstractPersistenceAlgorithm<Simplex> algo;
BarcodeCollection<Double> ints=null;
float[][] intervals;
int num_pts = 0;
PFont ft;
float int_max = 0;
float max = 0;
int toggle=1;

void settings() {
  fullScreen();
}

void setup() {
  ft = createFont("Courier", 16, true);
  textFont(ft, 14);
  float xs = width;
  float ys = height;
  background(255);
  line(xs/2, 0, xs/2, 0.8*ys);
  line(0, 0.8*ys, xs, 0.8*ys);
  fill(0);
  draw_instructions(0.08*xs, 0.8*ys, xs, ys);
  init_box(xs/2, 0, xs, 0.8*ys);
  fill(0);
  resetPoints();
  setupVRS();    
  pts1 = load_pts(1);
  pts2 = load_pts(0);
}

void draw() {
//  background(255);
  fill(255);
  rect(0, 0, width/2, 0.8*height);
  stroke(0);
  fill(0);
  
    for(Simplex s : vrs) {
    double fv = fc.getFiltrationValue(vrs.getFiltrationIndex(s));
    if(fv > f)
      continue;
    int[] ix;
    ix = s.getVertices();
    switch(s.getDimension()) { 
      case 0:
      if (toggle == 1){
            fill(255,255,255,0);
            ellipse((float)pts[ix[0]][0],(float)pts[ix[0]][1],(float)f,(float)f);
            fill(0);
            ellipse((float)pts[ix[0]][0],(float)pts[ix[0]][1],sizeX,sizeY);
            break;
      }
      if (toggle==-1){
        fill(0);
        ellipse((float)pts[ix[0]][0],(float)pts[ix[0]][1],sizeX,sizeY);
            break;
      }
      case 1:
        fill(0);
        line((float)pts[ix[0]][0],(float)pts[ix[0]][1],
            (float)pts[ix[1]][0],(float)pts[ix[1]][1]);
        break;
      case 2:
        fill(0,0,255,20);
        triangle((float)pts[ix[0]][0],(float)pts[ix[0]][1],
            (float)pts[ix[1]][0],(float)pts[ix[1]][1],
            (float)pts[ix[2]][0],(float)pts[ix[2]][1]);
        break;
      default:
        continue;
    }
  }
}

//*****************************************
// Compute a new VietorisRipsStream
//*****************************************

void setupVRS() {
  vrs = Plex4.createVietorisRipsStream(pts,2,maxeps,1000);
  fc = vrs.getConverter();
  ints=null;
}

//*****************************************
// Reset the points buffer
//*****************************************

void resetPoints() {
      pts=new double[0][2];
      dragX=0;
      dragY=0;
      offsetX=0;
      offsetY=0;
      sizeX=5;
      sizeY=5;
      f = 10;
      eps = 10;
      maxeps = 300;
}

//*****************************************
// Display instructions at bottom
//*****************************************

void draw_instructions(float xa, float ya, float xb, float yb) {
  int h = 14;
  text("INSTRUCTIONS", xa+50, h+ya);
  text("click        -- adds a point ", xa +10, 3*h + ya); 
  text("SHIFT-click  -- removes a point", xa +10, 4*h + ya); 
  text("1-4          -- loads example data sets", xa+10, 5*h + ya);
  text("T, Y, U, I   -- saves current data set to 5, 6, 7, or 8", xa+10, 6*h + ya);
  text("5-8          -- loads saved data set", xa+10, 7*h + ya);  
  text("RIGHT        -- step Vietoris-Rips complex forward", xa+10, 8*h+ya);
  text("LEFT         -- step Vietoris-Rips complex back", xa+10, 9*h+ya);
  text("C            -- clear points", xa + 10, 10*h+ya);
  text("Q            -- quit", xa+10, 11*h+ya);
  text("S            -- toggles Vietoris-Rips circles", xa+10,12*h+ya);
  text("InteractiveJPDwB. Luke Wolcott. 2016.", xa+850, 11*h+ya);
}

//*****************************************
// Initialize barcode box
//*****************************************

void init_box(float xa, float ya, float xb, float yb) {
  float cx = (xa + xb)/2;
  float cy = (ya + yb)/2;
  fill(255);
  stroke(0,0,204);
  strokeWeight(1.5);
  rect(cx - 250, cy - 250, 500, 500);
  stroke(5);
  strokeWeight(1);
}

//*****************************************
// On mousepress, if within data box then add a point. 
//*****************************************

void mousePressed() {
  if (keyPressed && keyCode == SHIFT){ 
    for (int i = 0; i < pts.length; i++){
      if (sq((float) pts[i][0]-mouseX)+sq((float) pts[i][1]-mouseY) < 25){  
        pts = remove_row(i);
        setupVRS();
        draw_barcode();
        break;
      }
    }
  }
  else{ 
    if ((mouseX < width/2) && (mouseY < 0.8*height)) {
      double[] pt = new double[2];

      translate(dragX,dragY);
      translate(offsetX,offsetY);

      pt[0] = mouseX;
      pt[1] = mouseY;
      
      println(pt[0]+","+pt[1]);
      pts = (double[][]) append(pts,pt);
      setupVRS();
      draw_barcode();
    }
  }
}

//*****************************************
// On keypress:
//
// Q      -- quit
// C      -- clear points
// B      -- run homology computation and plot barcode
// S      -- toggles circles 
// LEFT   -- step Vietoris-Rips complex back
// RIGHT  -- step Vietoris-Rips complex forward
// 1-4    -- loads pre-stored data sets
// 5-8    -- loads data sets saved by the user
// T,Y,U,I-- saves data set in 5-8, respectively
//*****************************************

void keyPressed() {
  float x_size = width/2;
  float y_size = 0.8*height;
  switch(key) {
    case 'q':
    case 'Q':
      exit();
      break;
      
    case 'c':
    case 'C':
      resetPoints();
      setupVRS();
      draw_barcode();
      break;
      
    case 's':
    case 'S':
       toggle = -toggle;
      
    case 'b':
    case 'B':
      draw_barcode();
      break;
    
    case '1':                              
      pts = pts1;
      setupVRS();
      draw_barcode();
      break;

    case '2':                              
      pts = pts2;
      setupVRS();
      draw_barcode();
      break;
      
    case 't':
    case 'T':
      pts5 = pts;
      break;
    
    case 'y':
    case 'Y':
      pts6 = pts;
      break;
    
    case 'u':
    case 'U':
      pts7 = pts;
      break;
    
    case 'i':
    case 'I':
      pts8 = pts;
      break;    
         
    case CODED:
      switch(keyCode) {
        case RIGHT:
          f += eps;
          println(f+": "+eps);
          if (f>max)
            f=max;
          draw_barcode();
          break;
        case LEFT:
          f -= eps;
          println(f+": "+eps);
          if(f<0)
            f=0;
          draw_barcode();
          break;
    }    
  }
}

//*****************************************
// Load pre-saved data into tables so it is ready to be read with 1-4 and 5-8.
//*****************************************

double[][] load_pts(int n) {
  double[][] ptsn;
  Table tablen;
  if (n==0){
    ptsn = new double[0][2];
    return ptsn;
  }
  else if (n==1)
    tablen = loadTable("seed1_data.csv", "header");
  else if (n==2)
    tablen = loadTable("seed2_cross.csv", "header");
  else if (n==3)
    tablen = loadTable("seed3_cross_hole_005.csv", "header");
  else 
    tablen = loadTable("seed4_circle.csv", "header");              
  ptsn = new double[tablen.getRowCount()][2];
  for (TableRow row : tablen.rows()){
    ptsn[row.getInt("point_id")][0] = row.getDouble("X_value")*(width/2)*0.8+(width/2)*0.05;
    ptsn[row.getInt("point_id")][1] = row.getDouble("Y_value")*(0.8*height)*0.8+(0.8*height)*0.05;
  }
  return ptsn;
}

//*****************************************
// Convert the output of a Javaplex persistence interval calculation into a tidy array.
//*****************************************

float[][] ints_to_intervals(String s){    
  String[] lines = splitTokens(s, " [,)\n\r");  

//  print how many lines, and then print each line
//  println("there are " + lines.length + " lines");
//  for (int i=0; i<lines.length; i++){
//    println(lines[i]);
//  }

  // Count how many different dimensions and points.
  int num_dims = 0;
  for (int i=0; i<lines.length; i++){
    if (lines[i].equals("Dimension:")){
      num_dims = num_dims + 1;
    }
  }      

  // print how many dimensions there are, and how many intervals there are
  println("there are " + num_dims + " different dimensions");
  num_pts = lines.length/2 - num_dims;
  println("there are " + num_pts + " different intervals");            

  // Build array of dimension, start, end.
  intervals = new float[num_pts][3];
  int dim = 0;
  int pt_number=-1;
  for (int k = 0; k<lines.length/2; k++){
    if (lines[2*k].equals("Dimension:")){
      dim = int(lines[2*k+1]);
      //println("dim:" + dim);
    }
    else{
      pt_number = pt_number + 1;
      //println(pt_number);
      intervals[pt_number][0] = dim;
      intervals[pt_number][1] = float(lines[2*k]);
      intervals[pt_number][2] = float(lines[2*k+1]);
    }  
  }  

  // print array of intervals
  for( int j=0; j<num_pts; j++){
    println(intervals[j][0], intervals[j][1], intervals[j][2]);
  } 
  
  return intervals;
}


//*****************************************
// Draw the barcode corresponding to a tidy array of intervals.
//*****************************************

void array_to_barcode(float[][] intervals){
    int nrow = intervals.length;
       
  // Look through table and figure out where the dimension changes.  
  int[] spots = {0};
  for (int i=1; i<nrow; i = i + 1){
    if (intervals[i][0] > intervals[i-1][0] ){
      spots = (int[]) append(spots,i);
    }
  }
    
  // Figure out what those dimensions actually are.  
  int[] dims = {int(intervals[0][0])};
  for (int i=1; i<spots.length; i=i+1){
    dims = (int[]) append(dims, int(intervals[spots[i]][0]));
  }
  spots = (int[]) append(spots,nrow);
  
  // Figure out horizontal scale.
  int_max = 0;
  for (int i=0; i<nrow; i=i+1){
    if (intervals[i][2] > int_max){
      int_max = intervals[i][2];
    }
  }
  
  println("Max interval value is " + int_max + ".");
  println("Infinite lines go all the way to end.");
  max = int_max*(1.1);        // Rescale so max isn't cut off. 

  // Convert infinity to max length.
  for (int i=0; i<nrow; i=i+1){
    if(Float.isNaN(intervals[i][2])){
      intervals[i][2] = max;
    }
  }     
  
  // Start drawing lines.
  float a = 0.75*width - 250;
  float b = 0.4*height - 250;
  float spaces = nrow - dims.length + 2 + 2 +4*(dims.length-1);
  float incr = 500/spaces;
  float y=b;      
  textSize(10);
  fill(0);
  for (int j=0; j<dims.length; j=j+1){
    y = y + 2*incr;
    text("dim " + dims[j], a-40, y);
    for (int k=spots[j]; k<spots[j+1]; k=k+1){    
      float start = intervals[k][1];
      float finish = intervals[k][2];
      line(a + 500*(start/max), y, a + 500*(finish/max), y);
      y = y + incr;
    }
    if (j < (dims.length-1)){
      y = y + 1*incr;
      stroke(0,0,204);
      strokeWeight(1.5);
      line(a, y, a+500, y);        // Draws a full line to separate dimensions
      stroke(0);
      strokeWeight(1);
    }
  }
  text(int(max) + " = max", a+490, b+520);
  text(0, a-2, b+520);
  stroke(0,0,204);
  strokeWeight(1.5);      
  line(a, b+500, a, b+510);
  line(a+500, b+500, a+500, b+510);
  line(a + (500/max)*((float) f), b + 500, a + (500/max)*((float) f), b + 510);
  stroke(0);
  strokeWeight(1);

}
  
void draw_barcode(){
      // compute intervals
      algo = Plex4.getDefaultSimplicialAlgorithm(2);
      ints = algo.computeIntervals(vrs);
      //println(ints);
      String s = ints.toString();
      
      // convert intervals into tidy array
      intervals = ints_to_intervals(s);
      
      //initialize barcode region   
      fill(255);
      rect(width/2, 0, width, 0.8*height);
      init_box(width/2, 0, width, 0.8*height);      
      
      // quit if there are no points in the array
      if (intervals.length == 0) {
        return;
      }
      
      // draw barcode from tidy array
      array_to_barcode(intervals);
}
  
//*****************************************
// Remove a row from pts.
//*****************************************
  
double[][] remove_row(int r){
  double[][] temp = new double[0][2];
  for (int i=0; i<pts.length; i++){
    if (i != r)
      temp = (double[][]) append(temp, pts[i]);
  }
  return temp;
}
  
  
  
  
  
  
  
  
  
  
