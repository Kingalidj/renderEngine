import processing.core.*; 
import processing.data.*; 
import processing.event.*; 
import processing.opengl.*; 

import java.util.Arrays; 

import java.util.HashMap; 
import java.util.ArrayList; 
import java.io.File; 
import java.io.BufferedReader; 
import java.io.PrintWriter; 
import java.io.InputStream; 
import java.io.OutputStream; 
import java.io.IOException; 

public class renderEngine extends PApplet {


canvas c;
PImage texture1;
mesh meshCube = new mesh();
float time = 0;
matrix matProj = new matrix(4, 4);
ArrayList < triangle > trisToRaster;

PVector camera = new PVector(0, 0, 0);
PVector lookDir;
float rotVal;

public void setup() {
  
  c = new canvas(width, height, 5);
  float aR = height / width;
  matProj = makeProjection(90, aR, 0.1f, 1000);
  texture1 = loadImage("mario.jpg");
  //meshCube.loadObj("hand");
}

public void draw() {
  //background(200, 0, 100);
  drawMesh(meshCube);
  c.background(200, 0, 100);
  c.show();
  //time += 0.03;
}

public void keyPressed() {
  if (keyCode == UP || keyCode == CONTROL) camera.y += 0.2f;
  if (keyCode == DOWN || key == ' ') camera.y -= 0.2f;
  if (keyCode == LEFT) camera.x += 0.2f;
  if (keyCode == RIGHT) camera.x -= 0.2f;

  if (key == 'A' || key == 'a') rotVal += 0.05f;
  if (key == 'D' || key == 'd') rotVal -= 0.05f;
  if (key == 'W' || key == 'w') {
    PVector f = PVector.mult(lookDir, 0.2f);
    camera.add(f);
  }
  if (key == 'S' || key == 's') {
    PVector f = PVector.mult(lookDir, 0.2f);
    camera.sub(f);
  }
}

public void drawMesh(mesh obj) {

  PVector lightDir = new PVector(0, -1, -1);

  trisToRaster = new ArrayList < triangle > ();

  matrix matRotZ = makeYRotation(time * 1), matRotX = makeXRotation(time * 0.5f);
  matrix matTrans = makeTranslation(0, 0, 8);

  matrix matWorld = new matrix(4, 4);
  matWorld.set(1);
  matWorld = multMatrix(matRotZ, matRotX);
  matWorld = multMatrix(matWorld, matTrans);

  PVector up = new PVector(0, 1, 0);
  PVector target = new PVector(0, 0, 1);
  matrix matCameraRot = makeYRotation(rotVal);
  lookDir = matCameraRot.mult(target);
  target = PVector.add(camera, lookDir);
  matrix matCamera = pointAtMatrix(camera, target, up);
  matrix matView = rotTransInverse(matCamera);

  for (triangle tri: obj.tris) {
    triangle triProj, triTrans, triView;

    triTrans = tri.applyMatrix(matWorld);
    triTrans.t = tri.t;

    PVector normal, line1, line2, cameraRay;

    line1 = PVector.sub(triTrans.p[1], triTrans.p[0]);
    line2 = PVector.sub(triTrans.p[2], triTrans.p[0]);
    normal = line1.cross(line2).normalize();

    cameraRay = PVector.sub(triTrans.p[0], camera);
    if (normal.dot(cameraRay) < 0) {

      lightDir.normalize();
      float dp = normal.dot(lightDir);
      float col = map(dp, 0, 1, 0, 255);
      triView = triTrans.applyMatrix(matView);
      triView.shade = color(col);
      triView.t = triTrans.t;

      triangle[] clips = new triangle[2];
      clips = clipTri(new PVector(0, 0, 0.1f), new PVector(0, 0, 1), triView);

      for (triangle clip: clips) {

        triProj = clip.applyMatrix(matProj);
        triProj.shade = clip.shade;
        triProj.translate(1, 1, 0);
        triProj.scale(0.5f * width, 0.5f * height);
        triProj.shade = clip.shade;
        triProj.t = clip.t;

        trisToRaster.add(triProj);
      }
    }
  }

  if (trisToRaster.size() > 0) {
    triangle[] trisToRasterArr = trisToRaster.toArray(new triangle[trisToRaster.size()]);
    Arrays.sort(trisToRasterArr);

    for (triangle triToRaster: trisToRasterArr) {

      ArrayList < triangle > triList = new ArrayList < triangle > ();
      triList.add(triToRaster);
      int newTriangles = 1;

      for (int i = 0; i < 4; i++) {
        triangle[] clipped = new triangle[0];
        while (newTriangles > 0) {
          triangle test = triList.get(0);
          triList.remove(0);
          newTriangles--;

          switch (i) {
            case 0:
              clipped = clipTri(new PVector(0, 0, 0), new PVector(0, 1, 0), test);
              break;
            case 1:
              clipped = clipTri(new PVector(0, height - 1, 0), new PVector(0, -1, 0), test);
              break;
            case 2:
              clipped = clipTri(new PVector(0, 0, 0), new PVector(1, 0, 0), test);
              break;
            case 3:
              clipped = clipTri(new PVector(width - 1, 0, 0), new PVector(-1, 0, 0), test);
              break;
          }

          for (triangle t: clipped) triList.add(t);
        }
        newTriangles = triList.size();
      }

      for (triangle t: triList) {
        texturedTri(round(t.p[0].x), round(t.p[0].y), t.t[0].x, t.t[0].y, round(t.p[1].x), round(t.p[1].y), t.t[1].x, t.t[1].y, round(t.p[2].x), round(t.p[2].y), t.t[2].x, t.t[2].y, texture1);
        t.show();
      }
    }
  }
}

class triangle implements Comparable < triangle > {
  PVector[] p = new PVector[3];
  PVector[] t = new PVector[3];
  int shade = 255;

  triangle() {
    p[0] = new PVector(0, 0);
    p[1] = new PVector(0, 0);
    p[2] = new PVector(0, 0);

    t[0] = new PVector(0, 0);
    t[1] = new PVector(0, 0);
    t[2] = new PVector(0, 0);
  }

  triangle(triangle tri) {
    p[0] = tri.p[0].copy();
    p[1] = tri.p[1].copy();
    p[2] = tri.p[2].copy();
  }

  triangle(PVector a, PVector b, PVector c) {
    p[0] = a;
    p[1] = b;
    p[2] = c;
  }

  triangle(PVector a, PVector b, PVector c, PVector d, PVector e, PVector f) {
    p[0] = a;
    p[1] = b;
    p[2] = c;

    t[0] = d;
    t[1] = e;
    t[2] = f;
  }

  public triangle applyMatrix(matrix m) {
    triangle triProjected = new triangle();
    for (int i = 0; i < 3; i++) triProjected.p[i] = m.mult(p[i]);
    return triProjected;
  }

  public int compareTo(triangle tri) {
    float aZ = (p[0].z + p[1].z + p[2].z) / 3;
    float otherAZ = (tri.p[0].z + tri.p[1].z + tri.p[2].z) / 3;
    if (aZ < otherAZ) return 1;
    if (aZ > otherAZ) return -1;
    return 0;
  }

  public void translate(float a, float b, float c) {
    for (PVector p: p) {
      p.add(a, b, c);
    }
  }

  public void scale(float a, float b) {
    for (PVector p: p) {
      p.x *= a;
      p.y *= b;
    }
  }

  public void show() {
    //stroke(shade);
    // stroke(200, 0, 100);
    //fill(shade);
    c.setFill(red(shade), green(shade), blue(shade));
    c.setStroke(255, 255, 255);
    c.edge = true;
    c.fill = false;
    //line(p[0].x, p[0].y, p[1].x, p[1].y);
    //line(p[1].x, p[1].y, p[2].x, p[2].y);
    //line(p[2].x, p[2].y, p[0].x, p[0].y);
    //triangle(p[0].x, p[0].y, p[1].x, p[1].y, p[2].x, p[2].y);
    c.triangle(p[0].x, p[0].y, p[1].x, p[1].y, p[2].x, p[2].y);
  }
}

public void texturedTri(int x1, int y1, float u1, float v1, int x2, int y2, float u2, float v2, int x3, int y3, float u3, float v3, PImage tex) {
  tex.loadPixels();
  if (y2 < y1) {
    int iTemp = y1;
    y1 = y2;
    y2 = iTemp;
    iTemp = x1;
    x1 = x2;
    x2 = iTemp;
    float fTemp = u1;
    u1 = u2;
    u2 = fTemp;
    fTemp = v1;
    v1 = v2;
    v2 = fTemp;
  }

  if (y3 < y1) {
    int iTemp = y1;
    y1 = y3;
    y3 = iTemp;
    iTemp = x1;
    x1 = x3;
    x3 = iTemp;
    float fTemp = u1;
    u1 = u3;
    u3 = fTemp;
    fTemp = v1;
    v1 = v3;
    v3 = fTemp;
  }

  if (y3 < y2) {
    int iTemp = y2;
    y2 = y3;
    y3 = iTemp;
    iTemp = x2;
    x2 = x3;
    x3 = iTemp;
    float fTemp = u2;
    u2 = u3;
    u3 = fTemp;
    fTemp = v2;
    v2 = v3;
    v3 = fTemp;
  }

  int dy1 = y2 - y1;
  int dx1 = x2 - x1;
  float du1 = u2 - u1;
  float dv1 = v2 - v1;

  int dy2 = y3 - y1;
  int dx2 = x3 - x1;
  float du2 = u3 - u1;
  float dv2 = v3 - v1;

  float texU, texV;

  float daxStep = 0, dbxStep = 0, du1Step = 0, dv1Step = 0, du2Step = 0, dv2Step = 0;

  if (dy1 != 0) daxStep = dx1 / (float) abs(dy1);
  if (dy2 != 0) dbxStep = dx2 / (float) abs(dy2);

  if (dy1 != 0) du1Step = du1 / (float) abs(dy1);
  if (dy2 != 0) dv1Step = dv1 / (float) abs(dy1);

  if (dy1 != 0) du2Step = du2 / (float) abs(dy2);
  if (dy2 != 0) dv2Step = dv2 / (float) abs(dy2);

  if (dy1 != 0) {
    for (int i = y1; i <= y2; i++) {
      int ax = PApplet.parseInt(x1 + (float)(i - y1) * daxStep);
      int bx = PApplet.parseInt(x1 + (float)(i - y1) * dbxStep);

      float su = u1 + (float)(i - y1) * du1Step;
      float sv = v1 + (float)(i - y1) * dv1Step;

      float eu = u1 + (float)(i - y1) * du2Step;
      float ev = v1 + (float)(i - y1) * dv2Step;

      if (ax > bx) {
        int iTemp = ax;
        ax = bx;
        bx = iTemp;
        float fTemp = su;
        su = eu;
        eu = fTemp;
        fTemp = sv;
        sv = ev;
        ev = fTemp;
      }

      texU = su;
      texV = sv;

      float tStep = 1 / ((float)(bx - ax));
      //--------------------------------------------------------------------------------------------------------
      float t = 0;

      for (int j = ax; j < bx; j++) {
        texU = (1 - t) * su + t * eu;
        texV = (1 - t) * sv + t * ev;

        int indx = PApplet.parseInt(texU + texV * tex.width);
        int col = tex.pixels[indx];
        c.setFill(red(col), green(col), blue(col));
        c.point(j, i);
        t += tStep;

      }
    }
  }

  dy1 = y3 - y2;
  dx1 = x3 - x2;
  dv1 = v3 - v2;
  du1 = u3 - u2;

  if (dy1 != 0) daxStep = dx1 / (float) abs(dy1);
  if (dy2 != 0) dbxStep = dx2 / (float) abs(dy2);

  if (dy1 != 0) du1Step = du1 / (float) abs(dy1);
  if (dy2 != 0) dv1Step = dv1 / (float) abs(dy1);

  if (dy1 != 0) {
    for (int i = y2; i <= y3; i++) {
      int ax = PApplet.parseInt(x2 + (float)(i - y2) * daxStep);
      int bx = PApplet.parseInt(x1 + (float)(i - y1) * dbxStep);

      float su = u2 + (float)(i - y2) * du1Step;
      float sv = v2 + (float)(i - y2) * dv1Step;

      float eu = u1 + (float)(i - y1) * du2Step;
      float ev = v1 + (float)(i - y1) * dv2Step;

      if (ax > bx) {
        int iTemp = ax;
        ax = bx;
        bx = iTemp;
        float fTemp = su;
        su = eu;
        eu = fTemp;
        fTemp = sv;
        sv = ev;
        ev = fTemp;
      }

      texU = su;
      texV = sv;

      float tStep = 1 / ((float)(bx - ax));
      //--------------------------------------------------------------------------------------------------------
      float t = 0;

      for (int j = ax; j < bx; j++) {
        texU = (1 - t) * su + t * eu;
        texV = (1 - t) * sv + t * ev;

        int indx = PApplet.parseInt(texU + texV * tex.width);
        int col = tex.pixels[indx];
        c.setFill(red(col), green(col), blue(col));
        c.point(j, i);
        t += tStep;

      }
    }
  }
}

class mesh {
  ArrayList < triangle > tris = new ArrayList < triangle > ();

  mesh() {
    tris.add(new triangle(new PVector(0, 0, 0), new PVector(0, 1, 0), new PVector(1, 1, 0), new PVector(0, 1), new PVector(0, 0), new PVector(1, 0)));
    tris.add(new triangle(new PVector(0, 0, 0), new PVector(1, 1, 0), new PVector(1, 0, 0), new PVector(0, 1), new PVector(0, 0), new PVector(1, 0)));
    tris.add(new triangle(new PVector(1, 0, 0), new PVector(1, 1, 0), new PVector(1, 1, 1), new PVector(0, 1), new PVector(0, 0), new PVector(1, 0)));
    tris.add(new triangle(new PVector(1, 0, 0), new PVector(1, 1, 1), new PVector(1, 0, 1), new PVector(0, 1), new PVector(0, 0), new PVector(1, 0)));
    tris.add(new triangle(new PVector(1, 0, 1), new PVector(1, 1, 1), new PVector(0, 1, 1), new PVector(0, 1), new PVector(0, 0), new PVector(1, 0)));
    tris.add(new triangle(new PVector(1, 0, 1), new PVector(0, 1, 1), new PVector(0, 0, 1), new PVector(0, 1), new PVector(0, 0), new PVector(1, 0)));
    tris.add(new triangle(new PVector(0, 0, 1), new PVector(0, 1, 1), new PVector(0, 1, 0), new PVector(0, 1), new PVector(0, 0), new PVector(1, 0)));
    tris.add(new triangle(new PVector(0, 0, 1), new PVector(0, 1, 0), new PVector(0, 0, 0), new PVector(0, 1), new PVector(0, 0), new PVector(1, 0)));
    tris.add(new triangle(new PVector(0, 1, 0), new PVector(0, 1, 1), new PVector(1, 1, 1), new PVector(0, 1), new PVector(0, 0), new PVector(1, 0)));
    tris.add(new triangle(new PVector(0, 1, 0), new PVector(1, 1, 1), new PVector(1, 1, 0), new PVector(0, 1), new PVector(0, 0), new PVector(1, 0)));
    tris.add(new triangle(new PVector(1, 0, 1), new PVector(0, 0, 1), new PVector(0, 0, 0), new PVector(0, 1), new PVector(0, 0), new PVector(1, 0)));
    tris.add(new triangle(new PVector(1, 0, 1), new PVector(0, 0, 0), new PVector(1, 0, 0), new PVector(0, 1), new PVector(0, 0), new PVector(1, 0)));
  }

  public void loadObj(String name) {
    name = name + ".obj";
    String[] obj = loadStrings(name);
    tris = new ArrayList < triangle > ();
    ArrayList < PVector > vertices = new ArrayList < PVector > ();

    for (int i = 0; i < obj.length; i++) {
      String line = obj[i];

      if (line.charAt(0) == 'v') {
        line = line.replace("v ", "");
        float[] numbers = PApplet.parseFloat(split(line, " "));
        PVector v = new PVector(numbers[0], numbers[1], numbers[2]);
        vertices.add(v);
      }

      if (line.charAt(0) == 'f') {
        line = line.replace("f ", "");
        int[] f = PApplet.parseInt(split(line, " "));
        tris.add(new triangle(vertices.get(f[0] - 1), vertices.get(f[1] - 1), vertices.get(f[2] - 1)));
      }
    }
  }
}

public matrix makeTranslation(float x, float y, float z) {
  matrix m = new matrix(4, 4);
  m.m[0][0] = 1;
  m.m[1][1] = 1;
  m.m[2][2] = 1;
  m.m[3][3] = 1;
  m.m[3][0] = x;
  m.m[3][1] = y;
  m.m[3][2] = z;
  return m;
}

public matrix makeProjection(float FOVDeg, float aR, float zNear, float zFar) {
  matrix m = new matrix(4, 4);
  float FOVRad = 1 / tan(FOVDeg * 0.5f / 180 * PI);
  m.m[0][0] = aR * FOVRad;
  m.m[1][1] = FOVRad;
  m.m[2][2] = zFar / (zFar - zNear);
  m.m[3][2] = (-zFar * zNear) / (zFar - zNear);
  m.m[2][3] = 1;
  return m;
}

public matrix makeYRotation(float angle) {
  matrix m = new matrix(4, 4);
  m.m[0][0] = cos(angle);
  m.m[0][2] = sin(angle);
  m.m[2][0] = -sin(angle);
  m.m[1][1] = 1;
  m.m[2][2] = cos(angle);
  m.m[3][3] = 1;
  return m;
}

public matrix makeXRotation(float angle) {
  matrix m = new matrix(4, 4);
  m.m[0][0] = 1;
  m.m[1][1] = cos(angle);
  m.m[1][2] = sin(angle);
  m.m[2][1] = -sin(angle);
  m.m[2][2] = cos(angle);
  m.m[3][3] = 1;
  return m;
}

public matrix makeZRotation(float angle) {
  matrix m = new matrix(4, 4);
  m.m[0][0] = cos(angle);
  m.m[0][1] = sin(angle);
  m.m[1][0] = -sin(angle);
  m.m[1][1] = cos(angle);
  m.m[2][2] = 1;
  m.m[3][3] = 1;
  return m;
}

public matrix multMatrix(matrix m1, matrix m2) {
  if (m1.cols == m2.rows) {
    matrix m3 = new matrix(m1.rows, m2.cols);
    for (int i = 0; i < m1.rows; i++) {
      for (int j = 0; j < m2.cols; j++) {
        float sum = 0;
        for (int k = 0; k < m1.cols; k++) {
          sum += m1.m[i][k] * m2.m[k][j];
        }
        m3.m[i][j] = sum;
      }
    }
    return m3;
  } else {
    println("multiplication not possible! number of cols or rows is incorrect.");
    return null;
  }
}

public matrix pointAtMatrix(PVector pos, PVector target, PVector up) {
  PVector nForward = PVector.sub(target, pos);
  nForward.normalize();

  PVector a = PVector.mult(nForward, up.dot(nForward));
  PVector nUp = PVector.sub(up, a).normalize();

  PVector nRight = nUp.cross(nForward);

  matrix m = new matrix(4, 4);
  m.m[0][0] = nRight.x;
  m.m[1][0] = nUp.x;
  m.m[2][0] = nForward.x;
  m.m[3][0] = pos.x;

  m.m[0][1] = nRight.y;
  m.m[1][1] = nUp.y;
  m.m[2][1] = nForward.y;
  m.m[3][1] = pos.y;

  m.m[0][2] = nRight.z;
  m.m[1][2] = nUp.z;
  m.m[2][2] = nForward.z;
  m.m[3][2] = pos.z;

  return m;
}

public matrix rotTransInverse(matrix m1) {
  matrix m2 = new matrix(4, 4);
  m2.m[0][0] = m1.m[0][0];
  m2.m[0][1] = m1.m[1][0];
  m2.m[0][2] = m1.m[2][0];
  m2.m[0][3] = 0;
  m2.m[1][0] = m1.m[0][1];
  m2.m[1][1] = m1.m[1][1];
  m2.m[1][2] = m1.m[2][1];
  m2.m[1][3] = 0;
  m2.m[2][0] = m1.m[0][2];
  m2.m[2][1] = m1.m[1][2];
  m2.m[2][2] = m1.m[2][2];
  m2.m[2][3] = 0;
  m2.m[3][0] = -(m1.m[3][0] * m2.m[0][0] + m1.m[3][1] * m2.m[1][0] + m1.m[3][2] * m2.m[2][0]);
  m2.m[3][1] = -(m1.m[3][0] * m2.m[0][1] + m1.m[3][1] * m2.m[1][1] + m1.m[3][2] * m2.m[2][1]);
  m2.m[3][2] = -(m1.m[3][0] * m2.m[0][2] + m1.m[3][1] * m2.m[1][2] + m1.m[3][2] * m2.m[2][2]);
  m2.m[3][3] = 1;

  return m2;
}

class matrix {
  float[][] m;
  int cols, rows;

  matrix(int x, int y) {
    rows = x;
    cols = y;
    m = new float[x][y];
  }

  public PVector mult(PVector i) {
    PVector o = new PVector(0, 0, 0);
    o.x = i.x * m[0][0] + i.y * m[1][0] + i.z * m[2][0] + m[3][0];
    o.y = i.x * m[0][1] + i.y * m[1][1] + i.z * m[2][1] + m[3][1];
    o.z = i.x * m[0][2] + i.y * m[1][2] + i.z * m[2][2] + m[3][2];
    float w = i.x * m[0][3] + i.y * m[1][3] + i.z * m[2][3] + m[3][3];

    if (w != 0) {
      o.x /= w;
      o.y /= w;
      o.z /= w;
    }
    return o;
  }


  public void set(float a) {
    for (int i = 0; i < rows; i++)
      for (int j = 0; j < cols; j++) m[i][j] = a;
  }
}

public PVector intersectPlane(PVector planeP, PVector planeN, PVector a, PVector b) {
  planeN.normalize();
  float planeD = -planeN.dot(planeP);
  float ad = a.dot(planeN);
  float bd = b.dot(planeN);
  float t = (-planeD - ad) / (bd - ad);
  PVector line = PVector.sub(b, a);
  PVector q = line.mult(t);
  return PVector.add(a, q);
}

public float getIntersectDist(PVector planeP, PVector planeN, PVector a, PVector b) {
  planeN.normalize();
  float planeD = -planeN.dot(planeP);
  float ad = a.dot(planeN);
  float bd = b.dot(planeN);
  float t = (-planeD - ad) / (bd - ad);
  return t;
}

public float _dist(PVector p, PVector planeN, PVector planeP) {
  //return (dist(planeN.x, planeN.y, planeN.z, p.x, p.y, p.z) - planeN.dot(planeP));
  return (planeN.x * p.x + planeN.y * p.y + planeN.z * p.z - planeN.dot(planeP));
}

public triangle[] clipTri(PVector planeP, PVector planeN, triangle inTri) {
  planeN.normalize();
  triangle[] out = new triangle[0];

  PVector[] inPoint = new PVector[3];
  PVector[] outPoint = new PVector[3];
  int inPointCount = 0, outPointCount = 0;

  PVector[] inTex = new PVector[3];
  PVector[] outTex = new PVector[3];
  int inTexCount = 0, outTexCount = 0;

  float d0 = _dist(inTri.p[0], planeN, planeP);
  float d1 = _dist(inTri.p[1], planeN, planeP);
  float d2 = _dist(inTri.p[2], planeN, planeP);

  if (d0 >= 0) {
    inPoint[inPointCount++] = inTri.p[0];
    inTex[inTexCount++] = inTri.t[0];
  } else {
    outPoint[outPointCount++] = inTri.p[0];
    outTex[outTexCount++] = inTri.t[0];
  }

  if (d1 >= 0) {
    inPoint[inPointCount++] = inTri.p[1];
    inTex[inTexCount++] = inTri.t[1];
  } else {
    outPoint[outPointCount++] = inTri.p[1];
    outTex[outTexCount++] = inTri.t[1];
  }

  if (d2 >= 0) {
    inPoint[inPointCount++] = inTri.p[2];
    inTex[inTexCount++] = inTri.t[2];
  } else {
    outPoint[outPointCount++] = inTri.p[2];
    outTex[outTexCount++] = inTri.t[2];
  }

  if (inPointCount == 0) {
    out = new triangle[0];
    return out;
  }

  if (inPointCount == 3) {
    out = new triangle[1];
    out[0] = inTri;

    return out;
  }

  if (inPointCount == 1 && outPointCount == 2) {
    out = new triangle[1];
    out[0] = new triangle();

    out[0].p[0] = inPoint[0];
    out[0].p[1] = intersectPlane(planeP, planeN, inPoint[0], outPoint[0]);
    out[0].p[2] = intersectPlane(planeP, planeN, inPoint[0], outPoint[1]);
    out[0].shade = inTri.shade;

    out[0].t[0] = inTex[0];

    float t = getIntersectDist(planeP, planeN, inPoint[0], outPoint[0]);
    out[0].t[1] = outTex[0].sub(inTex[0]).mult(t).add(inTex[0]);

    t = getIntersectDist(planeP, planeN, inPoint[0], outPoint[1]);
    out[0].t[2] = outTex[1].sub(inTex[0]).mult(t).add(inTex[0]);

    return out;
  }

  if (inPointCount == 2 && outPointCount == 1) {
    out = new triangle[2];
    out[0] = new triangle();
    out[1] = new triangle();

    out[0].p[0] = inPoint[0];
    out[0].p[1] = inPoint[1];
    out[0].p[2] = intersectPlane(planeP, planeN, inPoint[0], outPoint[0]);

    out[1].p[0] = inPoint[1];
    out[1].p[1] = out[0].p[2];
    out[1].p[2] = intersectPlane(planeP, planeN, inPoint[1], outPoint[0]);

    out[0].shade = inTri.shade;
    out[1].shade = inTri.shade;

    out[0].t[0] = inTex[0];
    out[0].t[1] = inTex[1];
    out[1].t[0] = inTex[1];
    out[1].t[1] = out[0].t[2].copy();

    float t = getIntersectDist(planeP, planeN, inPoint[0], outPoint[0]);
    out[0].t[2] = outTex[0].sub(inTex[0]).mult(t).add(inTex[0]);

    t = getIntersectDist(planeP, planeN, inPoint[1], outPoint[0]);
    out[1].t[2] = outTex[0].sub(inTex[1]).mult(t).add(inTex[0]);

    return out;
  }

  out = new triangle[0];
  return out;
}



class canvas {
  int[][] grid;
  boolean[][] update;
  int rows, cols;
  float size;
  float res;
  int col = color(255, 255, 255);
  int stroke = color(0, 0, 0);
  int background = color(0, 0, 0);
  boolean edge = true;
  boolean fill = true;
  ArrayList < PVector > changes = new ArrayList < PVector > ();

  canvas(float y, float x, float s) {
    rows = floor(x / s);
    cols = floor(y / s);
    size = s;
    res = 1 / s;
    grid = new int[cols][rows];
    update = new boolean[cols][rows];

    for (int i = 0; i < cols; i++)
      for (int j = 0; j < rows; j++) {
        update[i][j] = true;
        grid[i][j] = 10000;
      }
  }

  public void background(float r, float g, float b) {
    background = color(r, g, b);
    for (int i = 0; i < cols; i++)
      for (int j = 0; j < rows; j++)
        if (grid[i][j] != background && !update[i][j]) {
          changes.add(new PVector(i, j));
          grid[i][j] = background;
          update[i][j] = true;
        }
  }

  public void show() {
    noStroke();
    for (PVector p: changes) {
      fill(grid[round(p.x)][round(p.y)]);
      rect(p.x * size, p.y * size, size, size);
    }
    changes = new ArrayList < PVector > ();
    for (int i = 0; i < cols; i++)
      for (int j = 0; j < rows; j++) update[i][j] = false;
  }

  public void circle(float x1, float x2, float r) {
    PVector pos = new PVector(x1, x2);
    float rRes = res * r;
    float aRes = r * (pow(res, -0.03f)) * 5;
    if (fill) {
      for (int i = 0; i < rRes; i++) {
        for (float j = 0; j < aRes; j++) {
          PVector p = PVector.fromAngle(TWO_PI * j / aRes).normalize().mult(r);
          p.mult(PApplet.parseFloat(i) / rRes);
          p.add(pos).div(size);
          p.x = round(p.x);
          p.y = round(p.y);
          if (p.x < cols && p.x >= 0 && p.y < rows && p.y >= 0) {
            update[round(p.x)][round(p.y)] = true;
            if (grid[round(p.x)][round(p.y)] != col) {
              changes.add(p.copy());
              grid[round(p.x)][round(p.y)] = col;
            }
          }
        }
      }
    }
    if (edge) {
      for (int i = 0; i < aRes; i++) {
        PVector p = PVector.fromAngle(TWO_PI * i / aRes).normalize().mult(r);
        p.add(pos).div(size);
        p.x = round(p.x);
        p.y = round(p.y);
        if (p.x < cols && p.x >= 0 && p.y < rows && p.y >= 0) {
          update[round(p.x)][round(p.y)] = true;
          if (grid[round(p.x)][round(p.y)] != stroke) {
            changes.add(p.copy());
            grid[round(p.x)][round(p.y)] = stroke;
          }
        }
      }
    }
  }

  public void rectangle(float x1, float y1, float w, float h) {
    PVector pos = new PVector(x1 - w / 2, y1 - h / 2);
    PVector line1 = new PVector(w, 0);
    PVector line2 = new PVector(0, h);

    float line1Res = res * line1.mag();
    float line2Res = res * line2.mag();

    if (fill) {
      for (int i = 0; i <= line1Res; i++) {
        PVector p1 = PVector.mult(line1, PApplet.parseFloat(i) / line1Res);
        for (int j = 0; j <= line2Res; j++) {
          PVector p2 = PVector.mult(line2, PApplet.parseFloat(j) / line2Res);
          PVector p = PVector.add(p1, p2).add(pos);
          p.div(size);
          p.x = round(p.x);
          p.y = round(p.y);
          if (p.x < cols && p.x >= 0 && p.y < rows && p.y >= 0) {
            update[round(p.x)][round(p.y)] = true;
            if (grid[round(p.x)][round(p.y)] != col) {
              changes.add(p.copy());
              grid[round(p.x)][round(p.y)] = col;
            }
          }
        }
      }
    }
    if (edge) {
      line(pos.x, pos.y, pos.x + w, pos.y);
      line(pos.x, pos.y, pos.x, pos.y + h);
      line(pos.x + w, pos.y + h, pos.x, pos.y + h);
      line(pos.x + w, pos.y + h, pos.x + w, pos.y);
    }
  }

  public void triangle(float x1, float y1, float x2, float y2, float x3, float y3) {
    PVector line1 = new PVector(x2 - x1, y2 - y1);
    PVector line2 = new PVector(x3 - x1, y3 - y1);
    PVector a = new PVector(x1, y1);
    float line1Res = res * line1.mag() * 1.5f;
    float line2Res = res * line2.mag() * 1.5f;

    if (fill) {
      for (int i = 0; i <= line1Res; i++) {
        PVector p1 = PVector.mult(line1, PApplet.parseFloat(i) / line1Res);
        for (int j = 0; j <= line2Res * (1 - PApplet.parseFloat(i) / line1Res); j++) {
          PVector p2 = PVector.mult(line2, PApplet.parseFloat(j) / line2Res);
          PVector p = PVector.add(p1, p2).add(a);
          p.div(size);
          p.x = round(p.x);
          p.y = round(p.y);
          if (p.x < cols && p.x >= 0 && p.y < rows && p.y >= 0) {
            update[round(p.x)][round(p.y)] = true;
            if (grid[round(p.x)][round(p.y)] != col) {
              changes.add(p.copy());
              grid[round(p.x)][round(p.y)] = col;
            }
          }
        }
      }
    }
    if (edge) {
      line(x1, y1, x2, y2);
      line(x2, y2, x3, y3);
      line(x3, y3, x1, y1);
    }
  }


  public void line(float x1, float y1, float x2, float y2) {
    PVector a = new PVector(x1, y1);
    PVector b = new PVector(x2, y2);
    PVector line = PVector.sub(b, a);

    float lineRes = res * line.mag();

    for (int i = 0; i <= lineRes; i++) {
      PVector p = PVector.mult(line, PApplet.parseFloat(i) / lineRes).add(a);
      p.div(size);
      p.x = round(p.x);
      p.y = round(p.y);
      if (p.x < cols && p.x >= 0 && p.y < rows && p.y >= 0) {
        update[round(p.x)][round(p.y)] = true;
        if (grid[round(p.x)][round(p.y)] != stroke) {
          grid[round(p.x)][round(p.y)] = stroke;
          changes.add(p.copy());
        }
      }
    }
  }

  public void point(float x, float y) {
    PVector p = new PVector(x, y);
    p.div(size);
    if (p.x < cols && p.x >= 0 && p.y < rows && p.y >= 0) {
      update[round(p.x)][round(p.y)] = true;
      if (grid[round(p.x)][round(p.y)] != col) {
        grid[round(p.x)][round(p.y)] = col;
        changes.add(new PVector(round(p.x), round(p.y)));
      }
    }
  }



  public void setFill(float r, float g, float b) {
    col = color(r, g, b);
  }

  public void setStroke(float r, float g, float b) {
    stroke = color(r, g, b);
  }
}
  public void settings() {  size(1200, 1200); }
  static public void main(String[] passedArgs) {
    String[] appletArgs = new String[] { "--present", "--window-color=#666666", "--hide-stop", "renderEngine" };
    if (passedArgs != null) {
      PApplet.main(concat(appletArgs, passedArgs));
    } else {
      PApplet.main(appletArgs);
    }
  }
}
