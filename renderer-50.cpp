#include <glm/glm.hpp>
#include "sdw/DrawingWindow.h"
#include "sdw/CanvasPoint.h"
#include "sdw/CanvasTriangle.h"
#include "sdw/Colour.h"
#include "sdw/ModelTriangle.h"
#include "sdw/RayTriangleIntersection.h"
#include "sdw/TextureMap.h"
#include "sdw/Utils.h"
#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <limits>
#include <cmath>
#include <tuple>
#include <array>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <filesystem>
std::unordered_map<std::string,Colour> loadMtlColours(
  const std::string &mtlPath
);
std::vector<ModelTriangle> loadObjModelTriangles(
  const std::string &objPath, 
  float scale, 
  const std::unordered_map<std::string,Colour> &palette
); 
static void drawTexturedTriangleWithDepth(
  DrawingWindow &window, 
  const CanvasTriangle &tri, 
  const std::array<TexturePoint,3> &uv, 
  const TextureMap &tex, 
  std::vector<float> &depthBuffer
);
static void drawStrokedTriangle(
  DrawingWindow &window, 
  const CanvasTriangle &tri, 
  const Colour &colour
);
//just draw a line between two points:(
static void drawLine(
  DrawingWindow &window, 
  const CanvasPoint &from, 
  const CanvasPoint &to, 
  const Colour &colour
); 
//its about calculate the light intensity at a point
static float calculateProximityLighting(
  const glm::vec3 &point, 
  const glm::vec3 &lightPos, 
  float lightIntensity, 
  float ambientLight
); 
//to judge if a point is occluded by other triangles?
static bool isOccluded(
  const glm::vec3 &origin, 
  const glm::vec3 &lightPos, 
  const std::vector<ModelTriangle> &tris, 
  size_t excludeIdx
); 
static float calculateDiffuseLighting(
  const glm::vec3 &point,
  const glm::vec3 &normal, 
  const glm::vec3 &lightPos, 
  float lightIntensity, 
  float ambientLight, 
  const std::vector<ModelTriangle> &tris, 
  size_t triIndex, 
  bool enableShadow
); 
static RayTriangleIntersection getClosestValidIntersection(
  const glm::vec3 &cameraPos, 
  const glm::vec3 &rayDir, 
  const std::vector<ModelTriangle> &tris
);
//ARGB colors (which used on W2lab and MOCK mid-term exam :))
//but not in the midterm exam :(
static inline uint32_t packARGB(const Colour &c) { 
  return (255u << 24) | (static_cast<uint32_t>(c.red) << 16) | 
         (static_cast<uint32_t>(c.green) << 8) |               
         (static_cast<uint32_t>(c.blue));                      
}
static void savePPM(DrawingWindow &window, const std::string &filename) {
    int width  = window.width;
    int height = window.height;

    std::ofstream ofs(filename, std::ios::binary);//here is for windows(trans to /n),but idk if doesnt work on linux);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            uint32_t argb = window.getPixelColour(x, y);   // window.pixels[y*width + x] another choice but failed
            uint8_t r = (argb >> 16) & 0xFF;
            uint8_t g = (argb >> 8)  & 0xFF;
            uint8_t b = (argb)       & 0xFF;

            ofs.write(reinterpret_cast<char*>(&r), 1);
            ofs.write(reinterpret_cast<char*>(&g), 1);
            ofs.write(reinterpret_cast<char*>(&b), 1);
        }
    }
}
//week2 task 6 the barycentric coordinates of a point in a triangle
static inline glm::vec3 barycentric(const glm::vec2 &a, const glm::vec2 &b, const glm::vec2 &c, const glm::vec2 &p) { 
  float denominator = (b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y); 
  if (std::abs(denominator) < 1e-8f) {
    return glm::vec3(-1.0f);
  }     
  float weightA = ((b.y - c.y) * (p.x - c.x) + (c.x - b.x) * (p.y - c.y)) / denominator; 
  float weightB = ((c.y - a.y) * (p.x - c.x) + (a.x - c.x) * (p.y - c.y)) / denominator; 
  float weightC = 1.0f - weightA - weightB;                                              
  return glm::vec3(weightA, weightB, weightC);                                           
}
// this is the function of Utils.cpp but cant run in the future traingle lab task(still keep it)
// glm::vec3 convertToBarycentricCoordinates(glm::vec2 v0, glm::vec2 v1, glm::vec2 v2, glm::vec2 r)
// {
//     glm::vec2 e0 = v1 - v0;
//     glm::vec2 e1 = v2 - v0;
//     glm::vec2 e2 = r - v0;
//     float d00 = glm::dot(e0, e0);
//     float d01 = glm::dot(e0, e1);
//     float d11 = glm::dot(e1, e1);
//     float d20 = glm::dot(e2, e0);
//     float d21 = glm::dot(e2, e1);
//     float denominator = d00 * d11 - d01 * d01;
//     float u = (d11 * d20 - d01 * d21) / denominator;
//     float v = (d00 * d21 - d01 * d20) / denominator;
//     float w = 1.0f - u - v;
//     return glm::vec3(u,v,w);
// }
 static inline glm::mat3 lookatorientation(const glm::vec3 &cameraPos,
                                                  const glm::vec3 &target,
                                                  const glm::vec3 &worldUp = glm::vec3(0,1,0)) {
   
   glm::vec3 Back = glm::normalize(cameraPos - target);
   //Just for fault-tolerant processing:
 //avoid the situation that the calculated vector length is close to zero
   if (glm::length(Back) < 1e-6f) {
    Back = glm::vec3(0, 0, 1);
  }
   glm::vec3 Right = glm::normalize(glm::cross(worldUp, Back));//up x back = right yes!
   if (glm::length(Right) < 1e-6f) {
    Right = glm::normalize(glm::cross(glm::vec3(0,0,1), Back));
  }
   glm::vec3 Up = glm::normalize(glm::cross(Back, Right)
  );
return glm::mat3(Right, Up, Back);
 }
 static inline void lookAt(
  glm::mat3 &cameraOrientation,
  const glm::vec3 &cameraPos, 
  const glm::vec3 &target) {
   cameraOrientation = lookatorientation(cameraPos, target);
 }
//its about project a vertex onto the canvas point(week4 task7)
 static inline CanvasPoint vertextoCanvasPoint(
  const glm::vec3 &cameraPos, float focalLength,
  const glm::vec3 &vertex, int width, int height, float canvasScale,
  const glm::mat3 &cameraOrientation) {
    //glm::vec3 cameraToVertex = cameraPos - vertex; //is a previous wrong
     glm::vec3 cameraToVertex = vertex - cameraPos;                        
     glm::vec3 adjustedVector = cameraToVertex * cameraOrientation;    
     float x = adjustedVector.x;
     float y = adjustedVector.y;
     float z = adjustedVector.z;
     float u = -focalLength * (x / z);
   //The calculation of u includes a negative sign. This is because in camera space, z is negative (vertices lie in front), resulting in negative values for x/z.
   //Multiplying by -f yields a positive u. This ensures the X coordinate on the image plane aligns directionally with the X coordinate in world space.
     float v =  focalLength * (y / z);
     float cx = u * canvasScale + static_cast<float>(width) * 0.5f;
     float cy = v * canvasScale + static_cast<float>(height) * 0.5f;
     CanvasPoint p(cx, cy);
     //z-buffer
     p.depth = -1.0f / z;
     return p;
 }

//its a failed version code(but i  still to keep it for unexpected  situation)
// static inline CanvasPoint projectVertexOntoCanvasPoint(const glm::vec3 &cameraPos, float focalLength,
//                                                        const glm::vec3 &vertex, int width, int height, float canvasScale) {
//   float x = vertex.x - cameraPos.x;
//   float y = vertex.y - cameraPos.y;
//   float z = vertex.z - cameraPos.z;
//   float u = -focalLength * (x / z);
//   float v = focalLength * (y / z);
//   float cx = u * canvasScale + static_cast<float>(width) * 0.5f;
//   float cy = v * canvasScale + static_cast<float>(height) * 0.5f;
//   CanvasPoint p(cx, cy);
//   p.depth = -1.0f / z; 
//   return p;
// }

//used z-buffer to draw tri
static void filledtriangle(
  DrawingWindow &window, 
  const CanvasTriangle &tri,
  const Colour &colour, 
  std::vector<float> &depthBuffer) {
        int width = window.width;
        int height = window.height;
        //make a border of x and y (it will reduce the calculation haha)
        float minXfloat = std::min({tri[0].x, tri[1].x, tri[2].x});
        float maxXfloat = std::max({tri[0].x, tri[1].x, tri[2].x});
        float minYfloat = std::min({tri[0].y, tri[1].y, tri[2].y});
        float maxYfloat = std::max({tri[0].y, tri[1].y, tri[2].y});
        //change the float x and y to int x and y(for example let minxfloat=2.3f then minx=2...  )
        int minY = std::max(0, static_cast<int>(std::floor(minYfloat)));
        int maxY = std::min(height - 1, static_cast<int>(std::ceil(maxYfloat)));
        int minX = std::max(0, static_cast<int>(std::floor(minXfloat)));
        int maxX = std::min(width - 1, static_cast<int>(std::ceil(maxXfloat)));

        auto intersectEdge = [](
          const CanvasPoint &p0, 
          const CanvasPoint &p1, 
          float y) {
            float y0 = p0.y, y1 = p1.y;
            float ymin = std::min(y0, y1);
            float ymax = std::max(y0, y1);
            if (y <= ymin || y >= ymax){
              return std::pair<bool, std::pair<float, float>>(false, {0.0f, 0.0f});
  }
    //< and > results in redundant calculations of intersection points when scan lines pass precisely through vertices
    //(leading to superfluous pixels during filling).
            float t = (y - y0) / (y1 - y0);
            float x = p0.x + t * (p1.x - p0.x);
            float z = p0.depth + t * (p1.depth - p0.depth);
            return std::pair<bool, std::pair<float, float>>(true, {x, z});
  };

  for(int y = minY; y <= maxY; ++y){
    //float centerY = static_cast<float>(y);   no 0.5f added but seems like no change
    float centerY = static_cast<float>(y) + 0.5f;
    std::pair<bool, std::pair<float, float>> e01 = intersectEdge(tri[0], tri[1], centerY);
    std::pair<bool, std::pair<float, float>> e12 = intersectEdge(tri[1], tri[2], centerY);
    std::pair<bool, std::pair<float, float>> e20 = intersectEdge(tri[2], tri[0], centerY);
    std::vector<std::pair<float, float>> intersectionPoints;
    if(e01.first){
      intersectionPoints.push_back(e01.second);
    }
    if(e12.first){
      intersectionPoints.push_back(e12.second);
    }
    if(e20.first){
      intersectionPoints.push_back(e20.second);
      // if(intersectionPoints.size() == 3){
      //   break;
       }
      if(intersectionPoints.size() != 2){
        continue;
      }
      if(intersectionPoints[0].first > intersectionPoints[1].first){
        std::swap(intersectionPoints[0], intersectionPoints[1]);
      }
      float LeftX = intersectionPoints[0].first;
      //float RightX = intersectionPoints[1].first;
      float LeftZ = intersectionPoints[0].second;
      float RightX = intersectionPoints[1].first;
      float RightZ = intersectionPoints[1].second;
      
      if (LeftX == RightX){
        continue;
      }
      // int StartX = std::max(0, static_cast<int>(std::ceil(LeftX)));
      // int EndX = std::min(width - 1, static_cast<int>(std::floor(RightX)));(failed)
      int xStart = std::max(0, static_cast<int>(std::ceil(std::min(LeftX, RightX))));
      int xEnd = std::min(width - 1, static_cast<int>(std::floor(std::max(LeftX, RightX))));

       for (int x = xStart; x <= xEnd; ++x) {
      float xc = static_cast<float>(x) + 0.5f;
      float t = (xc - LeftX) / (RightX - LeftX);
      float invZ = LeftZ + t * (RightZ - LeftZ);
      int idx = y * width + x;
      if (invZ >= depthBuffer[idx]) {
        depthBuffer[idx] = invZ;
        window.setPixelColour(x, y, packARGB(colour));
      }
    

    }
  }
}
static void rotateCameraPositionAroundCentre(
  glm::vec3 &cameraPos, 
  const glm::vec3 &centre, 
  float angleRadians, 
  char axis) { 
    glm::vec3 distance = cameraPos - centre;
    float cosr = cos(angleRadians);
    float sinr = sin(angleRadians);
    glm::mat3 R(1.0f);
    if (axis == 'y') {
      R = glm::mat3(
      glm::vec3( cosr, 0, -sinr),
      glm::vec3( 0, 1,  0),
      glm::vec3( sinr, 0,  cosr)
    );
  } else if (axis == 'x') {
    R=glm::mat3(
      glm::vec3(1, 0, 0),
      glm::vec3(0,  cosr,  sinr),
      glm::vec3(0, -sinr,  cosr)
    );
  } else {//no z for the workbook
    return;
  }
  distance = R * distance;
  cameraPos = centre + distance;
}
static glm::mat3 XAxisRotationMatrix(float angleRadians) {
    float cosr = cos(angleRadians);
    float sinr = sin(angleRadians);
    return glm::mat3(
        glm::vec3(1.0f, 0.0f, 0.0f),
        glm::vec3(0.0f, cosr, sinr),
        glm::vec3(0.0f, -sinr, cosr)
    );
}

static glm::mat3 YAxisRotationMatrix(float angleRadians) {
    float c = cos(angleRadians);
    float s = sin(angleRadians);//c and s is easy than sinr and cosr write.......
    return glm::mat3(
        glm::vec3(c, 0.0f, -s),
        glm::vec3(0.0f, 1.0f, 0.0f),
        glm::vec3(s, 0.0f, c)
    );
}



static void rotateCameraOrientation(glm::mat3 &orientation, float angleRadians, char axis) {
    glm::mat3 rotationMatrix(1.0f);  

    if (axis == 'y') {
        
        rotationMatrix = YAxisRotationMatrix(angleRadians);
    } else if (axis == 'x') {
        
        rotationMatrix = XAxisRotationMatrix(angleRadians);
    } else {
        return;  //also no z for the workbook
    }

    
    orientation = rotationMatrix * orientation;
}
// draw a line between two points(Wireframe 3D scene rendering)
static void drawLine(DrawingWindow &window, const CanvasPoint &from, const CanvasPoint &to, const Colour &colour) {
    float xDiff = to.x - from.x;
    float yDiff = to.y - from.y;
    float numberOfSteps = std::max(std::abs(xDiff), std::abs(yDiff));
    if (numberOfSteps == 0) return;

    float xStepSize = xDiff / numberOfSteps;
    float yStepSize = yDiff / numberOfSteps;

    for (float i = 0.0f; i <= numberOfSteps; i++) {
        int x = static_cast<int>(std::round(from.x + (xStepSize * i)));
        int y = static_cast<int>(std::round(from.y + (yStepSize * i)));
        if (x >= 0 && x < window.width && y >= 0 && y < window.height) {
            window.setPixelColour(x, y, packARGB(colour));
        }
    }
}
// draw color of wireframe 3d
static void drawStrokedTriangle(DrawingWindow &window, const CanvasTriangle &tri, const Colour &colour) {
    drawLine(window, tri[0], tri[1], colour);
    drawLine(window, tri[1], tri[2], colour);
    drawLine(window, tri[2], tri[0], colour);
}
static float calculateProximityLighting(const glm::vec3 &point, const glm::vec3 &lightPos, float lightIntensity, float ambientLight) {
    float distance = glm::length(lightPos - point);
    if (distance < 0.01f) distance = 0.01f; // if r = 0 than will boom :(
// I = P / (4πr²)
    float brightness = lightIntensity / (4.0f * 3.14159265f * distance * distance);

    brightness = std::max(ambientLight, brightness);//just avoid too dark when is to far
    brightness = std::min(1.0f, brightness);
    brightness = std::max(0.0f, brightness);
return brightness;
}

// its to judge if a point is occluded by other triangles......
// origin: point position
// lightPos: light source position
// tris: all triangles
// excludeIdx: avoid self-occlusion(theres a bug when not exclude itself)
static bool isOccluded(const glm::vec3 &origin, const glm::vec3 &lightPos,
                       const std::vector<ModelTriangle> &tris, size_t excludeIdx) {
    const float eps = 1e-4f;
    glm::vec3 dir = glm::normalize(lightPos - origin);
    float maxDistance = glm::length(lightPos - origin);
    std::vector<ModelTriangle> filtered;
    filtered.reserve(tris.size());
    for (size_t i = 0; i < tris.size(); ++i) {
        if (i == excludeIdx) continue;
        filtered.push_back(tris[i]);
    }
    RayTriangleIntersection hit = getClosestValidIntersection(origin + eps * dir, dir, filtered);
    return hit.distanceFromCamera < maxDistance - eps;
}

// Diffuse lighting (proximity and angle-of-incidence)
// point: point position 
// normal: point normal
// lightIntensity: light source intensity
// ambientLight: minimum ambient light level
// triIndex: index of the triangle containing the point
static float calculateDiffuseLighting(const glm::vec3 &point, const glm::vec3 &normal,
                                      const glm::vec3 &lightPos, float lightIntensity,
                                      float ambientLight, const std::vector<ModelTriangle> &tris,
                                      size_t triIndex, bool enableShadow) {
    //same as calculateProximityLighting
    float distance = glm::length(lightPos - point);
    if (distance < 0.01f) distance = 0.01f;
    float distanceAttenuation = lightIntensity / (4.0f * 3.14159265f * distance * distance);

    // here is to calculate angle-of-incidence
    glm::vec3 lightDir = glm::normalize(lightPos - point);
    float angleOfIncidence = std::max(0.0f, glm::dot(glm::normalize(normal), lightDir));

    //Shadow calculation
    float shadowFactor = 1.0f;
    if (enableShadow && isOccluded(point, lightPos, tris, triIndex)) {
        shadowFactor = 0.0f; //completely in shadow
    }

    // final brightness calculation
    float brightness = std::max(ambientLight, distanceAttenuation * angleOfIncidence * shadowFactor);
    brightness = std::min(1.0f, brightness);
    brightness = std::max(0.0f, brightness);

    return brightness;
}
// soft shadow (tryfor 8 samples)
//also try for 16x but its too slow
static float calculateDiffuseLightingSoft(const glm::vec3 &point, const glm::vec3 &normal,
                                          const glm::vec3 &lightPos, float lightIntensity,
                                          float ambientLight, const std::vector<ModelTriangle> &tris,
                                          size_t triIndex) {
    glm::vec3 n = glm::normalize(normal);

    
    const int SAMPLE_COUNT = 8;      
    const float RADIUS = 0.12f;      

    float totalContribution = 0.0f;

    for (int i = 0; i < SAMPLE_COUNT; ++i) {
        //to make sure the total angle is 360
        float angle = (2.0f * 3.14159265f) * (static_cast<float>(i) / SAMPLE_COUNT);
        glm::vec3 samplePos = lightPos + RADIUS * glm::vec3(std::cos(angle), 0.0f, std::sin(angle));

        
        if (isOccluded(point, samplePos, tris, triIndex)) {
            continue;
        }

        glm::vec3 L = glm::normalize(samplePos - point);
        float ndotl = std::max(0.0f, glm::dot(n, L));

        float distance = glm::length(samplePos - point);
        if (distance < 0.01f) distance = 0.01f;
        float distanceAttenuation = lightIntensity / (4.0f * 3.14159265f * distance * distance);

        
        totalContribution += distanceAttenuation * ndotl;
        //to plus all thr light sources(8 light of 360)
    }

    // every light have same weight
    float avgContribution = 0.0f;
    if (SAMPLE_COUNT > 0) {
        avgContribution = totalContribution / static_cast<float>(SAMPLE_COUNT);
    }

    
    float brightness = std::max(ambientLight, avgContribution);
    brightness = std::min(1.0f, brightness);
    brightness = std::max(0.0f, brightness);
    return brightness;
}
//here is for specular lighting and phong(from w7)
static float gSpecularStrength = 1.0f;  //maybe smaller is better like 0.5f?  
static float gSpecularExponent = 256.0f;   
static float calculateSpecularLighting(
    const glm::vec3 &point,      
    const glm::vec3 &normal,     
    const glm::vec3 &lightPos,   
    const glm::vec3 &cameraPos,  
    float lightIntensity,        
    float specExponent           
) {
    glm::vec3 N = glm::normalize(normal);
    glm::vec3 L = glm::normalize(lightPos - point);
    glm::vec3 V = glm::normalize(cameraPos - point);
    glm::vec3 R = glm::reflect(-L, N);
    float rv = std::max(0.0f, glm::dot(R, V));
    float distance = glm::length(lightPos - point);
    if (distance < 0.01f) distance = 0.01f;
    float attenuation = lightIntensity / (4.0f * 3.14159265f * distance * distance);
    float spec = attenuation * std::pow(rv, specExponent);
    spec *= gSpecularStrength;
    return spec;
}
static RayTriangleIntersection getClosestValidIntersection(
  const glm::vec3 &cameraPos,
   const glm::vec3 &rayDir, 
   const std::vector<ModelTriangle> &tris) {
//here is three variables to store the best intersection info
    float bestT = std::numeric_limits<size_t>::max();//initially set to infinity
    size_t bestIdx = std::numeric_limits<size_t>::max();//initially set to max size_t
    RayTriangleIntersection best;//default constructor
    for (size_t i = 0; i < tris.size(); ++i) {
        const ModelTriangle &tri = tris[i];
        glm::vec3 e0 = tri.vertices[1] - tri.vertices[0];
        glm::vec3 e1 = tri.vertices[2] - tri.vertices[0];
        //E0 = P1 − P0,E1 = P2 − P0 for workbook
        glm::vec3 SPVector = cameraPos - tri.vertices[0];
        glm::mat3 DEMatrix(-rayDir, e0, e1);
        float det = glm::determinant(DEMatrix);
        if (std::abs(det) < 1e-8f) {
          continue;
        }
        glm::vec3 sol = glm::inverse(DEMatrix) * SPVector;
        float t = sol.x, u = sol.y, v = sol.z;
        //must to avoid some bug(like must be in front of the camera,must at 0-1,must in the triangle)
        if (t <= 0.0f) continue;
        if (u < 0.0f || u > 1.0f) continue;
        if (v < 0.0f || v > 1.0f) continue;
        if ((u + v) > 1.0f) continue;
        if (t < bestT) {
            glm::vec3 p = cameraPos + t * rayDir;
            best = RayTriangleIntersection(p, t, tri, i);
            bestT = t;
            bestIdx = i;
        }
    }
    if (bestIdx == std::numeric_limits<size_t>::max()) {
        return RayTriangleIntersection(glm::vec3(0.0f), std::numeric_limits<float>::max(), ModelTriangle(), std::numeric_limits<size_t>::max());
    }
  return best;
}
// maximum number of ray bounces for reflections(i think 4 is enough)
const int maxbounces = 4;
static Colour traceRay(
    const glm::vec3 &rayOrigin, 
    const glm::vec3 &rayDir,
    const std::vector<ModelTriangle> &tris, 
    const glm::vec3 &lightPos, 
    float lightIntensity, 
    float ambientLight,
    int bounceDepth,
    bool softShadowEnabled 
);
static glm::vec3 calculateRayDirection(int x, int y, int width, int height, 
                                       float focalLength, float canvasScale,
                                       const glm::mat3 &cameraOrientation) {
    float u = (static_cast<float>(x) + 0.5f - static_cast<float>(width) * 0.5f) / canvasScale;
    float v = (static_cast<float>(y) + 0.5f - static_cast<float>(height) * 0.5f) / canvasScale;
    // u = -focalLength * (cam.x / cam.z)
    // v =  focalLength * (cam.y / cam.z)
    // u = -focalLength * (cam.x / -focalLength) => cam.x = u
    // v =  focalLength * (cam.y / -focalLength) => cam.y = -v 
    //older version(not wrong i think but didnt work)
    glm::vec3 pixel_cam_space(u, -v, -focalLength); 
    glm::mat3 cameraToWorld = glm::transpose(cameraOrientation);
    glm::vec3 rayDir_world_space = glm::normalize(pixel_cam_space * cameraToWorld);
    return rayDir_world_space;
}
static Colour traceRay(
    const glm::vec3 &rayOrigin, 
    const glm::vec3 &rayDir,
    const std::vector<ModelTriangle> &tris, 
    const glm::vec3 &lightPos, 
    float lightIntensity, 
    float ambientLight,
    int bounceDepth,
    bool softShadowEnabled
) {
    if (bounceDepth > maxbounces) {
        return Colour(0, 0, 0); 
    }
    RayTriangleIntersection hit = getClosestValidIntersection(rayOrigin, rayDir, tris);
    if (!std::isfinite(hit.distanceFromCamera)) {
        return Colour(0, 0, 0); 
    }
    bool isMirror = (hit.intersectedTriangle.colour.name == "Yellow"); 
    if (isMirror) {
        //
        glm::vec3 Ri = rayDir; 
        glm::vec3 N = glm::normalize(hit.intersectedTriangle.normal);   
        if (glm::dot(Ri, N) > 0.0f) {
            N = -N;
        }
        glm::vec3 Rr = Ri - 2.0f * N * glm::dot(Ri, N); 
        //Rr​=Ri​−2(Ri​⋅N)N
        //////////it is important cuz of self-intersection
        glm::vec3 offsetOrigin = hit.intersectionPoint + (Rr * 1e-4f);
        return traceRay(offsetOrigin, Rr, tris, lightPos, lightIntensity, ambientLight, bounceDepth + 1, softShadowEnabled);
    } else {
        float brightness;
        if (softShadowEnabled) {
            brightness = calculateDiffuseLightingSoft(
                hit.intersectionPoint,
                hit.intersectedTriangle.normal,
                lightPos,
                lightIntensity,
                ambientLight,
                tris,
                hit.triangleIndex
            );
        } else {
            brightness = calculateDiffuseLighting(
                hit.intersectionPoint,
                hit.intersectedTriangle.normal,
                lightPos,
                lightIntensity,
                ambientLight,
                tris,
                hit.triangleIndex,
                true 
            );
        }      
        Colour baseColour = hit.intersectedTriangle.colour;
Colour litColour = baseColour;
        litColour.red = std::clamp(static_cast<int>(std::round(static_cast<float>(litColour.red) * brightness)), 0, 255);
        litColour.green = std::clamp(static_cast<int>(std::round(static_cast<float>(litColour.green) * brightness)), 0, 255);
        litColour.blue = std::clamp(static_cast<int>(std::round(static_cast<float>(litColour.blue) * brightness)), 0, 255);
        return litColour;
    }
}
enum ModelChoice { CornellModel, SphereModel, LogoModel};
static glm::vec3 computeSceneCentre(const std::vector<ModelTriangle> &tris) {
  float inf = std::numeric_limits<float>::max();
  glm::vec3 minv(inf, inf, inf);
  glm::vec3 maxv(-inf, -inf, -inf);
  //its to organize the scene and compute the centre point
  //its a stupid wrong like:
  //glm::vec3 minv(-inf, -inf, -inf);
 // glm::vec3 maxv(inf, inf, inf);
  for (const auto &mt : tris) {
    for (int i = 0; i < 3; ++i) {
      const glm::vec3 &v = mt.vertices[i];
      if (v.x < minv.x) minv.x = v.x;
      if (v.y < minv.y) minv.y = v.y;
      if (v.z < minv.z) minv.z = v.z;
      if (v.x > maxv.x) maxv.x = v.x;
      if (v.y > maxv.y) maxv.y = v.y;
      if (v.z > maxv.z) maxv.z = v.z;
    }
  }
  return (minv + maxv) * 0.5f;
  //same as previous version:
//centre.x = (min(x) + max(x)) / 2
//centre.y = (min(y) + max(y)) / 2
//centre.z = (min(z) + max(z)) / 2

}
///here is to load different obj model
static void loadModel(ModelChoice choice, std::vector<ModelTriangle> &tris, std::unordered_map<std::string, Colour> &palette) {
  if (choice == CornellModel) {
    palette = loadMtlColours("models/textured-cornell-box.mtl");
    tris = loadObjModelTriangles("models/textured-cornell-box.obj", 0.35f, palette);
  } else if(choice == SphereModel){
    palette.clear();
    tris = loadObjModelTriangles("resources/sphere.obj", 1.0f, palette);
  }else if(choice == LogoModel){
      palette.clear();
     tris = loadObjModelTriangles("models/logo.obj", 0.01f, palette);
   }
}
int main(int argc, char *argv[]) {
  const int width = 400;
  const int height = 300;
  DrawingWindow window(width, height, false);
  //camera parameters
  glm::vec3 cameraPos(0.0f, 0.0f, 4.0f);
  const float focalLength = 2.0f;
  const float canvasScale = 160.0f;

  //many main modes!!!!!!
  bool wireframeMode = false;
  bool proximityLightingMode = false;
  bool diffuseLightingMode = false;
  bool gouraudMode = false;
  bool textureEnabled = true;
  bool softShadowEnabled = false;
  bool rayTracingMode = false; 
  bool specularEnabled = false;  
  bool phongMode = false;  
   bool demoMode = false;     
  int  frameNumber = 0; 
  bool recording = true;
  // light parameters
  glm::vec3 lightPos(0.0f, 0.8f, 0.5f);  
  float lightIntensity = 20.0f;          
  float ambientLight = 0.2f;              

  std::unordered_map<std::string, Colour> palette;
  std::vector<ModelTriangle> tris;
  loadModel(CornellModel, tris, palette);
  TextureMap cobblesTex("models/texture.ppm");
  //TextureMap logoTex("models/texture2.ppm");//not work
  Colour cobblesFallback( 0, 255, 0);
  glm::vec3 sceneCentre = computeSceneCentre(tris);
  glm::mat3 cameraOrientation(1.0f);
  bool orbitEnabled = true;
  ModelChoice modelChoice = CornellModel;

  std::vector<float> depthBuffer(static_cast<size_t>(width * height), 0.0f);

//ok here is to judje whether the mode choose
  while (true) {
    window.clearPixels();
    std::fill(depthBuffer.begin(), depthBuffer.end(), 0.0f);
//     if (demoMode) {
//     if (frameNumber <= 200) {
//         wireframeMode = true;
//         diffuseLightingMode = false;
//         phongMode = false;
//     }
//     else if (frameNumber <= 400) {
//         wireframeMode = false;
//         diffuseLightingMode = true;
//         phongMode = false;
//     }
//     else if (frameNumber <= 600) {
//         wireframeMode = false;
//         diffuseLightingMode = false;
//         phongMode = true;
//     }
//     else {
        
//         frameNumber = 0;
//     }
// }//ok, here just a test i do for choose wireframe, diffuse, and phong mode,
if (demoMode) {// ok, this mode probly give up cuz i cnat find a perfact frame to show all the details
  //for each mode.......(for example we can see wireframe mode at the face but diffuselighting at the backside so it cant show the details well)
     
    if (frameNumber <= 100) {
        wireframeMode = true;
        proximityLightingMode = false;
        diffuseLightingMode = false;
        phongMode = false;
    }
    else if (frameNumber <= 200) {
        wireframeMode = false;
        proximityLightingMode = true;
        diffuseLightingMode = false;
        phongMode = false;
    }
    else if (frameNumber <= 300) {
        wireframeMode = false;
        diffuseLightingMode = false;
        phongMode = true;
    }
    else {
        
        frameNumber = 0;
    }
}   
    if (orbitEnabled) { 
      rotateCameraPositionAroundCentre(cameraPos, sceneCentre, 0.01f, 'y');
      lookAt(cameraOrientation, cameraPos, sceneCentre);
    }
    glm::mat3 activeOrientation = cameraOrientation;
    if (rayTracingMode) {
        
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                
                glm::vec3 rayDir = calculateRayDirection(x, y, width, height, focalLength, canvasScale, activeOrientation);
                
               
                Colour pixelColour = traceRay(
                    cameraPos, 
                    rayDir, 
                    tris, 
                    lightPos, 
                    lightIntensity, 
                    ambientLight, 
                    0, 
                    softShadowEnabled
                );
                
                
                window.setPixelColour(x, y, packARGB(pixelColour));
            }
        }

    } else {//Rasterization
        
        for (size_t triIndex = 0; triIndex < tris.size(); ++triIndex) {
           auto &mt = tris[triIndex];
          if (modelChoice == SphereModel) {
    mt.colour = Colour(255, 0, 0);
    mt.colour.name = "RedSphere";
}
          CanvasTriangle tri;
          for (int i = 0; i < 3; ++i) {
            tri[i] = vertextoCanvasPoint(cameraPos, focalLength, mt.vertices[i], width, height, canvasScale, activeOrientation);
          }

          //wireframe
          //else 
          if (wireframeMode) {
            //drawStrokedTriangle(window, tri, Colour(255,255,255));
            //here is whithe line of wireframe,but i think maybe mt.color is like more beautiful? :)
            drawStrokedTriangle(window, tri, mt.colour);
          } else if (proximityLightingMode) {//proximity lighting
            
            std::array<glm::vec3,3> wv{mt.vertices[0], mt.vertices[1], mt.vertices[2]};
            std::array<TexturePoint,3> uv{mt.texturePoints[0], mt.texturePoints[1], mt.texturePoints[2]};
            int width = window.width;
            int height = window.height;
            float minXf = std::min({tri[0].x, tri[1].x, tri[2].x});
            float maxXf = std::max({tri[0].x, tri[1].x, tri[2].x});
            float minYf = std::min({tri[0].y, tri[1].y, tri[2].y});
            float maxYf = std::max({tri[0].y, tri[1].y, tri[2].y});
            int minY = std::max(0, static_cast<int>(std::floor(minYf)));
            int maxY = std::min(height - 1, static_cast<int>(std::ceil(maxYf)));
            int minX = std::max(0, static_cast<int>(std::floor(minXf)));
            int maxX = std::min(width - 1, static_cast<int>(std::ceil(maxXf)));
            glm::vec2 a(tri[0].x, tri[0].y);
            glm::vec2 b(tri[1].x, tri[1].y);
            glm::vec2 c(tri[2].x, tri[2].y);
            for (int y = minY; y <= maxY; ++y) {
              float yc = static_cast<float>(y) + 0.5f;
              for (int x = minX; x <= maxX; ++x) {
                float xc = static_cast<float>(x) + 0.5f;
                glm::vec3 w = barycentric(a, b, c, glm::vec2(xc, yc));
                if (w.x < 0.0f || w.y < 0.0f || w.z < 0.0f) continue;
                float invZ = w.x * tri[0].depth + w.y * tri[1].depth + w.z * tri[2].depth;
                int idx = y * width + x;
                if (invZ > depthBuffer[idx]) {
                  float w0 = w.x * tri[0].depth;
                  float w1 = w.y * tri[1].depth;
                  float w2 = w.z * tri[2].depth;
                  float sumw = w0 + w1 + w2;
                  if (sumw <= 0.0f) continue;
                  glm::vec3 p = (wv[0] * w0 + wv[1] * w1 + wv[2] * w2) / sumw;
                  float brightness = calculateProximityLighting(p, lightPos, lightIntensity, ambientLight);
                  float uoz = w.x * (uv[0].x * tri[0].depth) + w.y * (uv[1].x * tri[1].depth) + w.z * (uv[2].x * tri[2].depth);
                  float voz = w.x * (uv[0].y * tri[0].depth) + w.y * (uv[1].y * tri[1].depth) + w.z * (uv[2].y * tri[2].depth);
                  float u = uoz / invZ;
                  float v = voz / invZ;
                  float uvSum = std::abs(uv[0].x) + std::abs(uv[0].y) + std::abs(uv[1].x) + std::abs(uv[1].y) + std::abs(uv[2].x) + std::abs(uv[2].y);
                  Colour baseColour = mt.colour;
                  if (!textureEnabled && uvSum > 1e-6f && mt.colour.name == std::string("Cobbles")) {
                    baseColour = cobblesFallback;
                  }
                  if (textureEnabled && uvSum > 1e-6f) {
                    int tx = std::clamp(static_cast<int>(u * static_cast<float>(cobblesTex.width - 1)), 0, static_cast<int>(cobblesTex.width - 1));
                    int ty = std::clamp(static_cast<int>((1.0f - v) * static_cast<float>(cobblesTex.height - 1)), 0, static_cast<int>(cobblesTex.height - 1));
                    uint32_t px = cobblesTex.pixels[ty * cobblesTex.width + tx];
                    baseColour.red = static_cast<uint8_t>((px >> 16) & 0xFF);
                    baseColour.green = static_cast<uint8_t>((px >> 8) & 0xFF);
                    baseColour.blue = static_cast<uint8_t>(px & 0xFF);
                  }
                  Colour litColour = baseColour;
                  litColour.red = std::clamp(static_cast<int>(std::round(static_cast<float>(litColour.red) * brightness)), 0, 255);
                  litColour.green = std::clamp(static_cast<int>(std::round(static_cast<float>(litColour.green) * brightness)), 0, 255);
                  litColour.blue = std::clamp(static_cast<int>(std::round(static_cast<float>(litColour.blue) * brightness)), 0, 255);
                  depthBuffer[idx] = invZ;
                  window.setPixelColour(x, y, packARGB(litColour));
                }
              }
            }
          } else if (diffuseLightingMode) {//diffuseLightingMode
            
            std::array<glm::vec3,3> wv{mt.vertices[0], mt.vertices[1], mt.vertices[2]};
            std::array<TexturePoint,3> uv{mt.texturePoints[0], mt.texturePoints[1], mt.texturePoints[2]};
            int width = window.width;
            int height = window.height;
            float minXf = std::min({tri[0].x, tri[1].x, tri[2].x});
            float maxXf = std::max({tri[0].x, tri[1].x, tri[2].x});
            float minYf = std::min({tri[0].y, tri[1].y, tri[2].y});
            float maxYf = std::max({tri[0].y, tri[1].y, tri[2].y});
            int minY = std::max(0, static_cast<int>(std::floor(minYf)));
            int maxY = std::min(height - 1, static_cast<int>(std::ceil(maxYf)));
            int minX = std::max(0, static_cast<int>(std::floor(minXf)));
            int maxX = std::min(width - 1, static_cast<int>(std::ceil(maxXf)));
            glm::vec2 a(tri[0].x, tri[0].y);
            glm::vec2 b(tri[1].x, tri[1].y);
            glm::vec2 c(tri[2].x, tri[2].y);
            for (int y = minY; y <= maxY; ++y) {
              float yc = static_cast<float>(y) + 0.5f;
              for (int x = minX; x <= maxX; ++x) {
                float xc = static_cast<float>(x) + 0.5f;
                glm::vec3 w = barycentric(a, b, c, glm::vec2(xc, yc));
                if (w.x < 0.0f || w.y < 0.0f || w.z < 0.0f) continue;
                float invZ = w.x * tri[0].depth + w.y * tri[1].depth + w.z * tri[2].depth;
                int idx = y * width + x;
                if (invZ >= depthBuffer[idx]) {
                  float w0 = w.x * tri[0].depth;
                  float w1 = w.y * tri[1].depth;
                  float w2 = w.z * tri[2].depth;
                  float sumw = w0 + w1 + w2;
                  if (sumw <= 0.0f) continue;
                  glm::vec3 p = (wv[0] * w0 + wv[1] * w1 + wv[2] * w2) / sumw;
                  float diffuseBrightness;
if (softShadowEnabled) { // soft shadow
    diffuseBrightness = calculateDiffuseLightingSoft(
        p, mt.normal,
        lightPos, lightIntensity, ambientLight,
        tris, triIndex
    );
} else { // hard shadow
    diffuseBrightness = calculateDiffuseLighting(
        p, mt.normal,
        lightPos, lightIntensity, ambientLight,
        tris, triIndex, true
    );
}
// its to judge if specular lighting is enabled
float specularBrightness = 0.0f;
if (specularEnabled) {
    specularBrightness = calculateSpecularLighting(
        p, //position(just for not forget)               
        mt.normal,       
        lightPos,
        cameraPos,
        lightIntensity,
        gSpecularExponent//256is good
    );
}
//diffuse + specular
float brightness = diffuseBrightness + specularBrightness;
brightness = std::min(1.0f, std::max(0.0f, brightness));
                  float uoz = w.x * (uv[0].x * tri[0].depth) + w.y * (uv[1].x * tri[1].depth) + w.z * (uv[2].x * tri[2].depth);
                  float voz = w.x * (uv[0].y * tri[0].depth) + w.y * (uv[1].y * tri[1].depth) + w.z * (uv[2].y * tri[2].depth);
                  float u = uoz / invZ;
                  float v = voz / invZ;
                  float uvSum = std::abs(uv[0].x) + std::abs(uv[0].y) + std::abs(uv[1].x) + std::abs(uv[1].y) + std::abs(uv[2].x) + std::abs(uv[2].y);
                  Colour baseColour = mt.colour;
                  if (!textureEnabled && uvSum > 1e-6f && mt.colour.name == std::string("Cobbles")) {
                    baseColour = cobblesFallback;
                  }
                  if (textureEnabled && uvSum > 1e-6f) {
                    int tx = std::clamp(static_cast<int>(u * static_cast<float>(cobblesTex.width - 1)), 0, static_cast<int>(cobblesTex.width - 1));
                    int ty = std::clamp(static_cast<int>((1.0f - v) * static_cast<float>(cobblesTex.height - 1)), 0, static_cast<int>(cobblesTex.height - 1));
                    uint32_t px = cobblesTex.pixels[ty * cobblesTex.width + tx];
                    baseColour.red = static_cast<uint8_t>((px >> 16) & 0xFF);
                    baseColour.green = static_cast<uint8_t>((px >> 8) & 0xFF);
                    baseColour.blue = static_cast<uint8_t>(px & 0xFF);
                  }
                  Colour litColour = baseColour;
                  litColour.red = std::clamp(static_cast<int>(std::round(static_cast<float>(litColour.red) * brightness)), 0, 255);
                  litColour.green = std::clamp(static_cast<int>(std::round(static_cast<float>(litColour.green) * brightness)), 0, 255);
                  litColour.blue = std::clamp(static_cast<int>(std::round(static_cast<float>(litColour.blue) * brightness)), 0, 255);
                  depthBuffer[idx] = invZ;
                  window.setPixelColour(x, y, packARGB(litColour));
                }
              }
            }
          } else if (gouraudMode) {//gouraud shading            
            std::array<glm::vec3,3> wv{mt.vertices[0], mt.vertices[1], mt.vertices[2]};
            std::array<TexturePoint,3> uv{mt.texturePoints[0], mt.texturePoints[1], mt.texturePoints[2]};
            int width = window.width;
            int height = window.height;
            float minXf = std::min({tri[0].x, tri[1].x, tri[2].x});
            float maxXf = std::max({tri[0].x, tri[1].x, tri[2].x});
            float minYf = std::min({tri[0].y, tri[1].y, tri[2].y});
            float maxYf = std::max({tri[0].y, tri[1].y, tri[2].y});
            int minY = std::max(0, static_cast<int>(std::floor(minYf)));
            int maxY = std::min(height - 1, static_cast<int>(std::ceil(maxYf)));
            int minX = std::max(0, static_cast<int>(std::floor(minXf)));
            int maxX = std::min(width - 1, static_cast<int>(std::ceil(maxXf)));
            glm::vec3 n0 = mt.normal;
            glm::vec3 n1 = mt.normal;
            glm::vec3 n2 = mt.normal;
            if (modelChoice == SphereModel) {
              glm::vec3 centre = sceneCentre;
              n0 = glm::normalize(wv[0] - centre);
              n1 = glm::normalize(wv[1] - centre);
              n2 = glm::normalize(wv[2] - centre);
            }
            float I0 = calculateDiffuseLighting(wv[0], n0, lightPos, lightIntensity, ambientLight, tris, triIndex, false);
            float I1 = calculateDiffuseLighting(wv[1], n1, lightPos, lightIntensity, ambientLight, tris, triIndex, false);
            float I2 = calculateDiffuseLighting(wv[2], n2, lightPos, lightIntensity, ambientLight, tris, triIndex, false);
            glm::vec2 a(tri[0].x, tri[0].y);
            glm::vec2 b(tri[1].x, tri[1].y);
            glm::vec2 c(tri[2].x, tri[2].y);
            for (int y = minY; y <= maxY; ++y) {
              float yc = static_cast<float>(y) + 0.5f;
              for (int x = minX; x <= maxX; ++x) {
                float xc = static_cast<float>(x) + 0.5f;
                glm::vec3 w = barycentric(a, b, c, glm::vec2(xc, yc));
                if (w.x < 0.0f || w.y < 0.0f || w.z < 0.0f) continue;
                float invZ = w.x * tri[0].depth + w.y * tri[1].depth + w.z * tri[2].depth;
                int idx = y * width + x;
                if (invZ >= depthBuffer[idx]) {
                  float w0 = w.x * tri[0].depth;
                  float w1 = w.y * tri[1].depth;
                  float w2 = w.z * tri[2].depth;
                  float sumw = w0 + w1 + w2;
                  if (sumw <= 0.0f) continue;
                  float brightness = (I0 * w0 + I1 * w1 + I2 * w2) / sumw;
                  float uoz = w.x * (uv[0].x * tri[0].depth) + w.y * (uv[1].x * tri[1].depth) + w.z * (uv[2].x * tri[2].depth);
                  float voz = w.x * (uv[0].y * tri[0].depth) + w.y * (uv[1].y * tri[1].depth) + w.z * (uv[2].y * tri[2].depth);
                  float u = uoz / invZ;
                  float v = voz / invZ;
                  float uvSum = std::abs(uv[0].x) + std::abs(uv[0].y) + std::abs(uv[1].x) + std::abs(uv[1].y) + std::abs(uv[2].x) + std::abs(uv[2].y);
                  Colour baseColour = mt.colour;
                  if (!textureEnabled && uvSum > 1e-6f && mt.colour.name == std::string("Cobbles")) {
                    baseColour = cobblesFallback;
                  }
                  if (textureEnabled && uvSum > 1e-6f) {
                    int tx = std::clamp(static_cast<int>(u * static_cast<float>(cobblesTex.width - 1)), 0, static_cast<int>(cobblesTex.width - 1));
                    int ty = std::clamp(static_cast<int>((1.0f - v) * static_cast<float>(cobblesTex.height - 1)), 0, static_cast<int>(cobblesTex.height - 1));
                    uint32_t px = cobblesTex.pixels[ty * cobblesTex.width + tx];
                    baseColour.red = static_cast<uint8_t>((px >> 16) & 0xFF);
                    baseColour.green = static_cast<uint8_t>((px >> 8) & 0xFF);
                    baseColour.blue = static_cast<uint8_t>(px & 0xFF);
                  }
                  Colour litColour = baseColour;
                  litColour.red = std::clamp(static_cast<int>(std::round(static_cast<float>(litColour.red) * brightness)), 0, 255);
                  litColour.green = std::clamp(static_cast<int>(std::round(static_cast<float>(litColour.green) * brightness)), 0, 255);
                  litColour.blue = std::clamp(static_cast<int>(std::round(static_cast<float>(litColour.blue) * brightness)), 0, 255);
                  depthBuffer[idx] = invZ;
                  window.setPixelColour(x, y, packARGB(litColour));
                }
              }
            }
          } 
          else if (phongMode) { //phong shading
            std::array<glm::vec3,3> wv{mt.vertices[0], mt.vertices[1], mt.vertices[2]};
            std::array<TexturePoint,3> uv{
        mt.texturePoints[0],
        mt.texturePoints[1],
        mt.texturePoints[2]
    };
            int width = window.width;
            int height = window.height;

            float minXf = std::min({tri[0].x, tri[1].x, tri[2].x});
            float maxXf = std::max({tri[0].x, tri[1].x, tri[2].x});
            float minYf = std::min({tri[0].y, tri[1].y, tri[2].y});
            float maxYf = std::max({tri[0].y, tri[1].y, tri[2].y});
            int minY = std::max(0, static_cast<int>(std::floor(minYf)));
            int maxY = std::min(height - 1, static_cast<int>(std::ceil(maxYf)));
            int minX = std::max(0, static_cast<int>(std::floor(minXf)));
            int maxX = std::min(width - 1, static_cast<int>(std::ceil(maxXf)));

            // ok,for this section i get two different normal for sphere and cornell box
            //yes its different between"Specular Lighting" and "phone shading"
            glm::vec3 n0 = mt.normal;
            glm::vec3 n1 = mt.normal;
            glm::vec3 n2 = mt.normal;
            if (modelChoice == SphereModel) {
              glm::vec3 centre = sceneCentre;
              n0 = glm::normalize(wv[0] - centre);
              n1 = glm::normalize(wv[1] - centre);
              n2 = glm::normalize(wv[2] - centre);
            }

            glm::vec2 a(tri[0].x, tri[0].y);
            glm::vec2 b(tri[1].x, tri[1].y);
            glm::vec2 c(tri[2].x, tri[2].y);

            for (int y = minY; y <= maxY; ++y) {
              float yc = static_cast<float>(y) + 0.5f;
              for (int x = minX; x <= maxX; ++x) {
                float xc = static_cast<float>(x) + 0.5f;
                glm::vec3 w = barycentric(a, b, c, glm::vec2(xc, yc));
                if (w.x < 0.0f || w.y < 0.0f || w.z < 0.0f) continue;

                float invZ = w.x * tri[0].depth + w.y * tri[1].depth + w.z * tri[2].depth;
                int idx = y * width + x;
                if (invZ < depthBuffer[idx]) continue;

    
                float w0 = w.x * tri[0].depth;
                float w1 = w.y * tri[1].depth;
                float w2 = w.z * tri[2].depth;
                float sumw = w0 + w1 + w2;
                if (sumw <= 0.0f) continue;

                glm::vec3 p = (wv[0] * w0 + wv[1] * w1 + wv[2] * w2) / sumw;              
                glm::vec3 n = glm::normalize((n0 * w0 + n1 * w1 + n2 * w2) / sumw);
                float uoz = w0 * uv[0].x + w1 * uv[1].x + w2 * uv[2].x;
                float voz = w0 * uv[0].y + w1 * uv[1].y + w2 * uv[2].y;
                float u = 0.0f, v = 0.0f;

              bool hasUV = false;
              
float uvSum = std::abs(uv[0].x) + std::abs(uv[0].y)
            + std::abs(uv[1].x) + std::abs(uv[1].y)
            + std::abs(uv[2].x) + std::abs(uv[2].y);
if (uvSum > 1e-6f && invZ > 0.0f) {
    u = uoz / invZ;
    v = voz / invZ;
    hasUV = true;
}

                //here is the phong lighting calculation
                glm::vec3 L = glm::normalize(lightPos - p);
                float ndotl = std::max(0.0f, glm::dot(n, L));

                float distance = glm::length(lightPos - p);
                if (distance < 0.01f) distance = 0.01f;
                float attenuation = lightIntensity / (4.0f * 3.14159265f * distance * distance);

                float diffuseTerm = ndotl * attenuation;

              
                float specularTerm = calculateSpecularLighting(
                    p, n, lightPos, cameraPos,
                    lightIntensity, gSpecularExponent
                );

                
                float brightness = std::min(1.0f, std::max(ambientLight, diffuseTerm));
                

                //here i choose red cuz white is so hard to see the effect
                Colour baseColour = mt.colour;
                if (textureEnabled && hasUV && mt.colour.name == std::string("Cobbles")) {
   int tx = std::clamp(
    static_cast<int>(u * (cobblesTex.width - 1)),
    0,
    static_cast<int>(cobblesTex.width - 1)
);
    int ty = std::clamp(
    static_cast<int>((1.0f - v) * (cobblesTex.height - 1)),
    0,
    static_cast<int>(cobblesTex.height - 1)
);

    uint32_t px = cobblesTex.pixels[ty * cobblesTex.width + tx];
    baseColour.red   = (px >> 16) & 0xFF;
    baseColour.green = (px >>  8) & 0xFF;
    baseColour.blue  =  px        & 0xFF;
}


if (!textureEnabled && hasUV && mt.colour.name == std::string("Cobbles")) {
    baseColour = cobblesFallback;
}
                if (modelChoice == SphereModel) {
                  baseColour = Colour(220, 40, 40);   
                }

                int r = std::clamp(
                    static_cast<int>(baseColour.red   * brightness + 255.0f * specularTerm), 0, 255);
                int g = std::clamp(
                    static_cast<int>(baseColour.green * brightness + 255.0f * specularTerm), 0, 255);
                int b = std::clamp(
                    static_cast<int>(baseColour.blue  * brightness + 255.0f * specularTerm), 0, 255);

                depthBuffer[idx] = invZ;
                window.setPixelColour(x, y, packARGB(Colour(r, g, b)));
              }
            }

          }else {//
            
            std::array<TexturePoint,3> uv{mt.texturePoints[0], mt.texturePoints[1], mt.texturePoints[2]};
            float uvSum = std::abs(uv[0].x) + std::abs(uv[0].y) + std::abs(uv[1].x) + std::abs(uv[1].y) + std::abs(uv[2].x) + std::abs(uv[2].y);
            if (textureEnabled && uvSum > 1e-6f) {
              drawTexturedTriangleWithDepth(window, tri, uv, cobblesTex, depthBuffer);
            } else {
              Colour drawColour = mt.colour;
              if (!textureEnabled && uvSum > 1e-6f && mt.colour.name == std::string("Cobbles")) {
                drawColour = cobblesFallback;
              }
              filledtriangle(window, tri, drawColour, depthBuffer);
            }
          }
        }
    }
    


    cameraOrientation = activeOrientation;

    
    SDL_Event event;
    while (window.pollForInputEvents(event)) {
      if (event.type == SDL_KEYDOWN) {
        switch (event.key.keysym.sym) {
          
          case SDLK_0: 
            rayTracingMode = true;
            wireframeMode = false;
            proximityLightingMode = false;
            diffuseLightingMode = false;
            gouraudMode = false;
            std::cout << "switch to ray" << std::endl;
            break;
          case SDLK_1:
            rayTracingMode = false; 
            wireframeMode = true;
            proximityLightingMode = false;
            diffuseLightingMode = false;
            gouraudMode = false;
            phongMode = false;
            std::cout << "switch to wire" << std::endl;
            break;
          case SDLK_2:
            rayTracingMode = false; 
            wireframeMode = false;
            proximityLightingMode = false;
            diffuseLightingMode = false;
            gouraudMode = false;
            phongMode = false;
            std::cout << "switch to rasterization" << std::endl;
            break;
          case SDLK_3:
            rayTracingMode = false; 
            wireframeMode = false;
            proximityLightingMode = true;
            diffuseLightingMode = false;
            gouraudMode = false;
            phongMode = false;
            std::cout << "switch to proximity light" << std::endl;
            break;
          case SDLK_4:
            rayTracingMode = false; 
            wireframeMode = false;
            proximityLightingMode = false;
            diffuseLightingMode = true;
            gouraudMode = false;
            phongMode = false;  
            std::cout << "switch to diffuse light" << std::endl;
            break;
          case SDLK_5:
            rayTracingMode = false; 
            wireframeMode = false;
            proximityLightingMode = false;
            diffuseLightingMode = false;
            gouraudMode = true;
            phongMode = false;  
            std::cout << "switch to gouraud" << std::endl;
            
            break;
          case SDLK_6:
            textureEnabled = !textureEnabled;
          
            break;
          case SDLK_e:               
            rayTracingMode = false;
            wireframeMode = false;
            proximityLightingMode = false;
            diffuseLightingMode = false;
            gouraudMode = false;
            phongMode = true;
            std::cout << "switch to Phong shading" << std::endl;
            break;
          
          case SDLK_w: cameraPos.z -= 0.1f; break;
          case SDLK_s: cameraPos.z += 0.1f; break;
          case SDLK_LEFT:  cameraPos.x -= 0.1f; break;
          case SDLK_RIGHT: cameraPos.x += 0.1f; break;
          case SDLK_UP:    cameraPos.y -= 0.1f; break;
          case SDLK_DOWN:  cameraPos.y += 0.1f; break;

          //a,d,e,f is special, youcan delete the lookAt function to have a free camera(use ijkl)
          case SDLK_a: rotateCameraPositionAroundCentre(cameraPos, sceneCentre, 0.05f, 'y');lookAt(cameraOrientation, cameraPos, sceneCentre); break;
          
          case SDLK_d: rotateCameraPositionAroundCentre(cameraPos, sceneCentre, -0.05f, 'y');lookAt(cameraOrientation, cameraPos, sceneCentre); break;
          case SDLK_r: rotateCameraPositionAroundCentre(cameraPos, sceneCentre, 0.05f, 'x'); lookAt(cameraOrientation, cameraPos, sceneCentre);break;
          case SDLK_f: rotateCameraPositionAroundCentre(cameraPos, sceneCentre, -0.05f, 'x'); lookAt(cameraOrientation, cameraPos, sceneCentre);break;

          
          case SDLK_j: rotateCameraOrientation(cameraOrientation, 0.05f, 'y'); break;
          case SDLK_l: rotateCameraOrientation(cameraOrientation, -0.05f, 'y'); break;
          case SDLK_i: rotateCameraOrientation(cameraOrientation, 0.05f, 'x'); break;
          case SDLK_k: rotateCameraOrientation(cameraOrientation, -0.05f, 'x'); break;

          
          case SDLK_o:
            orbitEnabled = !orbitEnabled;
            std::cout << "orbit: " << (orbitEnabled ? "on" : "off") << std::endl;
            break;

          // LookAt功能：相机朝向场景中心
          case SDLK_p:
            lookAt(cameraOrientation, cameraPos, sceneCentre);
            std::cout << "lookat" << std::endl;
            break;
          case SDLK_7:
            modelChoice = CornellModel;
            loadModel(modelChoice, tris, palette);
            sceneCentre = computeSceneCentre(tris);
            lightPos = glm::vec3(0.0f, 0.8f, 0.5f);//why i add the function cuz when i first select ball then box,the light still at ball position lol
lookAt(cameraOrientation, cameraPos, sceneCentre);
            std::cout << "Cornell Box" << std::endl;
            break;
          case SDLK_8:
            modelChoice = SphereModel;
            loadModel(modelChoice, tris, palette);
            sceneCentre = computeSceneCentre(tris);
            lightPos = glm::vec3(1.2f, 1.2f, 2.0f);
            lookAt(cameraOrientation, cameraPos, sceneCentre);
            std::cout << " ball" << std::endl;
            break;
            // case SDLK_z:
            // modelChoice = LogoModel;
            // loadModel(modelChoice, tris, palette);
            // sceneCentre = computeSceneCentre(tris);
            // std::cout << " ball" << std::endl;
            // break;

            //explain:at the begining i want to make a logo model which is different from cornell box and sphere,like7,8,z to switch
            //but later i find that its better to add logo into cornell box instead of making a separate model
            //so i comment the code above and implement the function below
            case SDLK_z: {

            std::cout << "Adding logo into Cornell Box…" << std::endl;

    
              std::unordered_map<std::string, Colour> logoPalette;
              std::vector<ModelTriangle> logoMesh =
                  loadObjModelTriangles("models/logo.obj", 0.002f, logoPalette);  

    
              for (auto &tri : logoMesh) {
                  for (int i = 0; i < 3; i++) {
                      tri.vertices[i] += glm::vec3(-1.0f, -1.0f, 0.5f);   //its to hard to find the position,maybe i need to add a controller lol
                  }
            }

    
    tris.insert(tris.end(), logoMesh.begin(), logoMesh.end());

    std::cout << "Logo added!" << std::endl;
    break;
}

          case SDLK_9:
            softShadowEnabled = !softShadowEnabled;
            std::cout << "shadow " << (softShadowEnabled ? "Soft Shadow (multi-point light)" : "Hard Shadow (single-point light)") << std::endl;
            break;
            case SDLK_v:
            recording = !recording;
            std::cout << "Recording: " << (recording ? "ON" : "OFF") << std::endl;
            break;
            

          
          case SDLK_t: lightPos.y += 0.1f;  break;
          case SDLK_g: lightPos.y -= 0.1f;   break;
          case SDLK_h: lightPos.x -= 0.1f;  break;
          case SDLK_y: lightPos.x += 0.1f;  break;
          case SDLK_u: lightPos.z -= 0.1f;  break;
          case SDLK_m: lightPos.z += 0.1f;  break;

          default: break;
        }
      }
    }

    window.renderFrame();
    //here is the recording part
    if (recording) {
        std::filesystem::create_directory("frames");
       std::ostringstream oss;
        oss << "frames/frame-" << std::setw(4) << std::setfill('0') << frameNumber << ".ppm";//here i use 0001-1000,maybe 00001-10000 later?
        std::string filename = oss.str();

        savePPM(window, filename);
        std::cout << "Saved frame: " << filename << std::endl;

        ++frameNumber;
    }

    SDL_Delay(16);
  }

  return 0;
}

//i think the most difficult part:some of the obj r defined like v/t/n or v//n or v/t
static std::pair<int,int> parseVertexTexIndexToken(const std::string &token) {
    auto parts = split(token, '/');
    if (parts.empty() || parts[0].empty()) throw std::invalid_argument("empty index");
    int vi = std::stoi(parts[0]);
    int ti = -1;
    if (parts.size() >= 2 && !parts[1].empty()) ti = std::stoi(parts[1]);
    return {vi, ti};
}

std::unordered_map<std::string, Colour> loadMtlColours(const std::string &mtlPath) {
    std::ifstream mtl(mtlPath);
    //privious code to judge if mtl is open
    // if (!mtl.is_open()) {
    //     throw std::runtime_error("Failed to open MTL: " + mtlPath);
    // }
    std::unordered_map<std::string, Colour> palette;
    std::string line, current;
    int lineNo = 0;
    while (std::getline(mtl, line)) {
        ++lineNo;
        if (line.empty()) continue;
        size_t start = line.find_first_not_of(" \t\r");
        if (start == std::string::npos) continue;
        std::string trimmed = line.substr(start);
        auto tokens = split(trimmed, ' ');
        tokens.erase(std::remove_if(tokens.begin(), tokens.end(), [](const std::string &s){ return s.empty(); }), tokens.end());
        for (auto &t : tokens) {
            if (!t.empty() && (t.back() == '\r' || t.back() == '\n')) t.pop_back();
        }
        if (tokens.empty()) continue;
        if (tokens[0] == "newmtl") {
            if (tokens.size() >= 2) current = tokens[1];
        } else if (tokens[0] == "Kd") {
            if (current.empty()) continue;
            if (tokens.size() >= 4) {
                float r = std::stof(tokens[1]);
                float g = std::stof(tokens[2]);
                float b = std::stof(tokens[3]);
                int ri = std::clamp(static_cast<int>(std::round(r * 255.0f)), 0, 255);
                int gi = std::clamp(static_cast<int>(std::round(g * 255.0f)), 0, 255);
                int bi = std::clamp(static_cast<int>(std::round(b * 255.0f)), 0, 255);
                Colour c(ri, gi, bi);
                c.name = current;
                palette[current] = c;
            }
        }
    }
    return palette;
}

std::vector<ModelTriangle> loadObjModelTriangles(const std::string &objPath, float scale, const std::unordered_map<std::string, Colour> &palette) {
    std::ifstream in(objPath);
    //judge if obj is open......
    // if (!in.is_open()) {
    //     throw std::runtime_error("Failed to open OBJ: " + objPath);
    // }

    std::vector<glm::vec3> vertices;
    std::vector<TexturePoint> texcoords;
    std::vector<ModelTriangle> triangles;
    Colour currentColour(255,255,255);
    currentColour.name = "DefaultWhite";

    std::string line; int lineNo = 0;
    while (std::getline(in, line)) {
        ++lineNo; if (line.empty()) continue;
        size_t start = line.find_first_not_of(" \t\r");
        if (start == std::string::npos) continue;
        std::string trimmed = line.substr(start);
        auto tokens = split(trimmed, ' ');
        tokens.erase(std::remove_if(tokens.begin(), tokens.end(), [](const std::string &s){ return s.empty(); }), tokens.end());
        for (auto &t : tokens) { if (!t.empty() && (t.back() == '\r' || t.back() == '\n')) t.pop_back(); }
        if (tokens.empty()) continue;
        const std::string &type = tokens[0];

        //at the beginning it error all time so i add a judge rules here,they can give me an output

        if (type == "v") {
            if (tokens.size() < 4) {// std::cerr << "Warning: line " << lineNo << " invalid vertex, skipped" << std::endl;
               continue; }
            float x = std::stof(tokens[1]); float y = std::stof(tokens[2]); float z = std::stof(tokens[3]);
            vertices.emplace_back(x * scale, y * scale, z * scale);
        } else if (type == "vt") {/////////////////////////
            if (tokens.size() < 3) { //std::cerr << "Warning: line " << lineNo << " invalid texcoord, skipped" << std::endl; 
              continue; }
            float u = std::stof(tokens[1]); float v = std::stof(tokens[2]);
            texcoords.emplace_back(u, v);
        } else if (type == "usemtl") {//////////////////////////////////
            if (tokens.size() >= 2) {
                auto it = palette.find(tokens[1]);
                if (it != palette.end()) currentColour = it->second; else {
                    //std::cerr << "Warning: material '" << tokens[1] << "' not in palette, using white" << std::endl;
                    currentColour = Colour(255,255,255); currentColour.name = "DefaultWhite";
                }
            }
        } else if (type == "f") {
            if (tokens.size() < 4) { //std::cerr << "Warning: line " << lineNo << " invalid face, skipped" << std::endl;
               continue; }
            std::vector<int> vi; vi.reserve(tokens.size()-1);
            std::vector<int> ti; ti.reserve(tokens.size()-1);
            for (size_t t = 1; t < tokens.size(); ++t) {//////////////////////
                try {
                    auto pr = parseVertexTexIndexToken(tokens[t]);
                    int vIndex = pr.first - 1;
                    int tIndex = (pr.second >= 1) ? pr.second - 1 : -1;
                    if (vIndex < 0 || static_cast<size_t>(vIndex) >= vertices.size()) throw std::out_of_range("face vertex index out of range");
                    vi.push_back(vIndex);
                    ti.push_back(tIndex);
                } catch (const std::exception &e) {/////////////////////////
                    //std::cerr << "Warning: line " << lineNo << " face index parse failed: " << e.what() << ", skipped" << std::endl;
                    vi.clear(); ti.clear(); break;
                }
            }
            if (vi.size() < 3) continue;
            for (size_t k = 2; k < vi.size(); ++k) {
                ModelTriangle tri;
                tri.vertices[0] = vertices[vi[0]];
                tri.vertices[1] = vertices[vi[k-1]];
                tri.vertices[2] = vertices[vi[k]];
                if (ti[0] >= 0 && ti[k-1] >= 0 && ti[k] >= 0 &&
                    static_cast<size_t>(ti[0]) < texcoords.size() &&
                    static_cast<size_t>(ti[k-1]) < texcoords.size() &&
                    static_cast<size_t>(ti[k]) < texcoords.size()) {
                    tri.texturePoints[0] = texcoords[ti[0]];
                    tri.texturePoints[1] = texcoords[ti[k-1]];
                    tri.texturePoints[2] = texcoords[ti[k]];
                }//////////////////////////
                tri.colour = currentColour;
                tri.normal = glm::normalize(glm::cross(tri.vertices[1] - tri.vertices[0], tri.vertices[2] - tri.vertices[0]));
                triangles.push_back(tri);
            }
        }
    }

    
    return triangles;
}

static void drawTexturedTriangleWithDepth(DrawingWindow &window, const CanvasTriangle &tri,/////////////////
                                          const std::array<TexturePoint,3> &uv, const TextureMap &tex,
                                          std::vector<float> &depthBuffer) {
  int width = window.width;
  int height = window.height;
  float minYf = std::min({tri[0].y, tri[1].y, tri[2].y});
  float maxYf = std::max({tri[0].y, tri[1].y, tri[2].y});
  int minY = std::max(0, static_cast<int>(std::floor(minYf)));
  int maxY = std::min(height - 1, static_cast<int>(std::ceil(maxYf)));

  auto interpEdge = [&](const CanvasPoint &p0, const CanvasPoint &p1, const TexturePoint &t0, const TexturePoint &t1, float y) {
    float y0 = p0.y, y1 = p1.y;
    float ymin = std::min(y0, y1);
    float ymax = std::max(y0, y1);
    if (y <= ymin || y >= ymax) return std::pair<bool, std::tuple<float,float,float,float>>(false, std::make_tuple(0,0,0,0));
    float tt = (y - y0) / (y1 - y0);
    float x = p0.x + tt * (p1.x - p0.x);
    float invZ = p0.depth + tt * (p1.depth - p0.depth);
    float uoz0 = t0.x * p0.depth; float uoz1 = t1.x * p1.depth;
    float voz0 = t0.y * p0.depth; float voz1 = t1.y * p1.depth;
    float uoz = uoz0 + tt * (uoz1 - uoz0);
    float voz = voz0 + tt * (voz1 - voz0);
    return std::pair<bool, std::tuple<float,float,float,float>>(true, std::make_tuple(x, invZ, uoz, voz));
  };

  for (int y = minY; y <= maxY; ++y) {
    float yc = static_cast<float>(y) + 0.5f;
    auto e01 = interpEdge(tri[0], tri[1], uv[0], uv[1], yc);
    auto e12 = interpEdge(tri[1], tri[2], uv[1], uv[2], yc);
    auto e20 = interpEdge(tri[2], tri[0], uv[2], uv[0], yc);
    std::vector<std::tuple<float,float,float,float>> xs;
    if (e01.first) xs.push_back(e01.second);
    if (e12.first) xs.push_back(e12.second);
    if (e20.first) xs.push_back(e20.second);
    if (xs.size() != 2) continue;
    if (std::get<0>(xs[0]) > std::get<0>(xs[1])) std::swap(xs[0], xs[1]);
    float xL = std::get<0>(xs[0]);
    float zL = std::get<1>(xs[0]);
    float uozL = std::get<2>(xs[0]);
    float vozL = std::get<3>(xs[0]);
    float xR = std::get<0>(xs[1]);
    float zR = std::get<1>(xs[1]);
    float uozR = std::get<2>(xs[1]);
    float vozR = std::get<3>(xs[1]);
    if (xR == xL) continue;
    int xStart = std::max(0, static_cast<int>(std::ceil(std::min(xL, xR))));
    int xEnd = std::min(width - 1, static_cast<int>(std::floor(std::max(xL, xR))));
    for (int x = xStart; x <= xEnd; ++x) {
      float xc = static_cast<float>(x) + 0.5f;
      float t = (xc - xL) / (xR - xL);
      float invZ = zL + t * (zR - zL);
      float uoz = uozL + t * (uozR - uozL);
      float voz = vozL + t * (vozR - vozL);
      int idx = y * width + x;
      if (invZ >= depthBuffer[idx]) {
        depthBuffer[idx] = invZ;
        float u = uoz / invZ;
        float v = voz / invZ;
        int tx = std::clamp(static_cast<int>(u * static_cast<float>(tex.width - 1)), 0, static_cast<int>(tex.width - 1));
        int ty = std::clamp(static_cast<int>((1.0f - v) * static_cast<float>(tex.height - 1)), 0, static_cast<int>(tex.height - 1));
        window.setPixelColour(x, y, tex.pixels[ty * tex.width + tx]);
      }
    }
  }
}