
using System;
using System.Buffers;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography.X509Certificates;
using Unity.VisualScripting;
using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.SocialPlatforms.GameCenter;
using UnityEngine.TextCore;
using UnityEngine.UIElements;

public class Terraformer : MonoBehaviour
{
    private MeshFilter meshFilter;
    private Mesh mesh;

    public int radius;

    private int xRadius, zRadius;
    public float noiseScale = 2f;

    public float height = 12f;


    public bool noise = true;

    public float groundHeight = 5f;

    [Range(0, 1)]
    public float groundSmoothing = 0.2f;

    int i = 0;

    public float torusRadius = 4.0f;

    public float waterFallsSeed = 1.23f;

    public float riverWidth = 1.0f;
    public float centerRiverCircleRadius = 1.0f;


    public float terrainNoiseScale = 0.25f;

    public MeshFilter water;


    // Start is called before the first frame update
    void Start()
    {
        meshFilter = GetComponent<MeshFilter>();
        mesh = new();

        // water.transform.localScale = new(radius/5, radius/5, radius/5);

        
        xRadius = radius;
        zRadius = radius;

        RenderMesh();
        RenderWaterMesh();

        meshFilter.mesh = mesh;
        lastSeed = waterFallsSeed;
    }


    private float lastSeed;
    // Update is called once per frame
    void Update()
    {
        if (waterFallsSeed != lastSeed) {
            RenderMesh();
            RenderWaterMesh();
        }

        lastSeed = waterFallsSeed;
    }

    void RenderMesh() {
        List<Vector3> verts = new();
        List<int> tris = new();

        Dictionary<Vector3, int> verts_to_index = new();

        List<int> boundaryVerts = new();

        List<int> topVerts = new();
        List<int> bottomVerts = new();
        List<int> leftVerts = new();
        List<int> rightVerts = new();

        for(int i = -xRadius + 1; i < xRadius; i++) {
            for(int j = -zRadius + 1; j < zRadius; j++) {


                Vector3[] newVerts = {
                    new(i - 0.5f, 0, j - 0.5f),
                    new(i + 0.5f, 0, j - 0.5f),
                    new(i + 0.5f, 0, j + 0.5f),
                    new(i - 0.5f, 0, j + 0.5f),
                };

                int[] indicies = {
                    0, 0, 0, 0
                };

                int ii = 0;
                foreach (Vector3 vertex in newVerts) {
                    int index;
                    if (verts_to_index.ContainsKey(vertex)) {
                        index = verts_to_index[vertex];
                        indicies[ii] = index;
                    } else {
                        index = verts.Count;
                        verts.Add(vertex);

                        verts_to_index.Add(vertex, index);
                        indicies[ii] = index;
                    }

                    if (vertex.x == xRadius - 0.5f && !rightVerts.Contains(index)) {
                        rightVerts.Add(index);
                    } else if (vertex.x == -xRadius + 0.5f && !leftVerts.Contains(index)) {
                        leftVerts.Add(index);
                    } else if (vertex.z == zRadius - 0.5f && !topVerts.Contains(index)) {
                        topVerts.Add(index);
                    } else if (vertex.z == -zRadius + 0.5f && !bottomVerts.Contains(index )) {
                        bottomVerts.Add(index);
                    }

                    ii++;
                }

                //create triangles
                tris.AddRange(new List<int> {
                    indicies[0], indicies[1], indicies[3],
                    indicies[1], indicies[2], indicies[3], 
                });
            }
        }


        for(int i = 0; i < verts.Count; i++) {
            // x or z value is this percentage to the edge
            float xRatioToEdge = Math.Abs(xRadius / verts[i].x);
            float zRatioToEdge = Math.Abs(zRadius / verts[i].z);

            float ratioToEdge = Mathf.Min(xRatioToEdge, zRatioToEdge);

            float length = (verts[i] * ratioToEdge).magnitude;

            float boundaryRatio = radius / length;

            // Debug.DrawLine(new Vector3(0, 0, 0), verts[i] * boundaryRatio);

            verts[i] *= boundaryRatio;
        }



        for(int i = 0; i < verts.Count; i++) {
            float xPercent = (verts[i].x / xRadius) / 2 + 0.5f;
            float zPercent = (verts[i].z / zRadius) / 2 + 0.5f;

            verts[i] += new Vector3(0, -EvaluateY(xPercent, zPercent) * height, 0);
        }


        rightVerts.Reverse();
        bottomVerts.Reverse();

        boundaryVerts.AddRange(topVerts);
        boundaryVerts.AddRange(rightVerts);
        boundaryVerts.AddRange(bottomVerts);
        boundaryVerts.AddRange(leftVerts);

        for(int i = 0; i < boundaryVerts.Count - 1; i++) {
            if(verts[boundaryVerts[i]] == verts[boundaryVerts[i+1]]) {
                 boundaryVerts.RemoveAt(i+1);
            }
        }

        if(verts[boundaryVerts[boundaryVerts.Count - 1]] == verts[boundaryVerts[0]]) {
            boundaryVerts.RemoveAt(boundaryVerts.Count - 1);
        }

        int offset = verts.Count;

        //first extrusion

        for(int i = 0; i < boundaryVerts.Count; i++) {
            verts.Add(verts[boundaryVerts[i]] + new Vector3(0f, groundHeight, 0f) + verts[boundaryVerts[i]].normalized * torusRadius);
        }
        for(int i = 0; i < boundaryVerts.Count - 1; i++) {
            tris.AddRange(new List<int> {
                    boundaryVerts[i], boundaryVerts[i+1], offset + i,
                    boundaryVerts[i + 1], offset + i + 1, offset + i,
                }
            );
        }
        tris.AddRange(new List<int> {
                    boundaryVerts[boundaryVerts.Count - 1], boundaryVerts[0], offset + boundaryVerts.Count - 1,
                    boundaryVerts[0], offset, offset + boundaryVerts.Count - 1,
            }
        );

        int secondExtrusionOffset = verts.Count;

        //second extrusion
        for(int i = 0; i < boundaryVerts.Count; i++) {
            verts.Add(verts[offset + i] + new Vector3(0, groundHeight, 0)  + verts[boundaryVerts[i]].normalized * torusRadius/4);
        }
        for(int i = 0; i < boundaryVerts.Count - 1; i++) {
            tris.AddRange(new List<int> {
                    offset + i, offset + i+1, secondExtrusionOffset + i,
                    offset +i + 1, secondExtrusionOffset + i + 1, secondExtrusionOffset + i,
                }
            );
        }
        tris.AddRange(new List<int> {
                    offset + boundaryVerts.Count - 1, offset, secondExtrusionOffset + boundaryVerts.Count - 1,
                    offset, secondExtrusionOffset, secondExtrusionOffset + boundaryVerts.Count - 1,
            }
        );










        verts_to_index.Clear();

        List<int> newBoundaryVerts = new();

        rightVerts.Clear();
        leftVerts.Clear();
        topVerts.Clear();
        bottomVerts.Clear();

        offset = verts.Count;

        List<int> terrainVerts = new();

        for(int i = -xRadius + 1; i < xRadius; i++) {
            for(int j = -zRadius + 1; j < zRadius; j++) {


                Vector3[] newVerts = {
                    new(i - 0.5f, 0, j - 0.5f),
                    new(i + 0.5f, 0, j - 0.5f),
                    new(i + 0.5f, 0, j + 0.5f),
                    new(i - 0.5f, 0, j + 0.5f),
                };

                int[] indicies = {
                    0, 0, 0, 0
                };

                int ii = 0;
                foreach (Vector3 vertex in newVerts) {
                    int index;
                    if (verts_to_index.ContainsKey(vertex)) {
                        index = verts_to_index[vertex];
                        indicies[ii] = index;
                    } else {
                        index = verts.Count;
                        verts.Add(vertex + new Vector3(0, 2*groundHeight, 0) + vertex.normalized*torusRadius);

                        verts_to_index.Add(vertex, index);
                        indicies[ii] = index;
                    }

                    if (vertex.x == xRadius - 0.5f && !rightVerts.Contains(index)) {
                        rightVerts.Add(index);
                    } else if (vertex.x == -xRadius + 0.5f && !leftVerts.Contains(index)) {
                        leftVerts.Add(index);
                    } else if (vertex.z == zRadius - 0.5f && !topVerts.Contains(index)) {
                        topVerts.Add(index);
                    } else if (vertex.z == -zRadius + 0.5f && !bottomVerts.Contains(index )) {
                        bottomVerts.Add(index);
                    } else {
                        terrainVerts.Add(index);
                    }

                    ii++;
                }

                //create triangles
                tris.AddRange(new List<int> {
                    indicies[0],  indicies[3], indicies[1],
                    indicies[1],  indicies[3], indicies[2],
                });
            }
        }



        rightVerts.Reverse();
        bottomVerts.Reverse();

        newBoundaryVerts.AddRange(topVerts);
        newBoundaryVerts.AddRange(rightVerts);
        newBoundaryVerts.AddRange(bottomVerts);
        newBoundaryVerts.AddRange(leftVerts);


        for(int i = offset; i < verts.Count; i++) {
            // x or z value is this percentage to the edge
            float xRatioToEdge = Math.Abs(xRadius / verts[i].x);
            float zRatioToEdge = Math.Abs(zRadius / verts[i].z);

            float ratioToEdge = Mathf.Min(xRatioToEdge, zRatioToEdge);

            float length = (verts[i] * ratioToEdge).magnitude;

            float boundaryRatio = radius / length;

            // Debug.DrawLine(new Vector3(0, 0, 0), verts[i] * boundaryRatio);

            verts[i] *= boundaryRatio;
        }

        for(int i = 0; i < boundaryVerts.Count - 1; i++) {
            tris.AddRange(new List<int> {
                secondExtrusionOffset + i,newBoundaryVerts[i + 1], newBoundaryVerts[i], 
                secondExtrusionOffset + i, secondExtrusionOffset + i + 1,  newBoundaryVerts[i + 1],
            });
            // Debug.DrawLine(verts[secondExtrusionOffset + i], verts[newBoundaryVerts[i]]);
        }
        tris.AddRange(new List<int> {
                secondExtrusionOffset + boundaryVerts.Count - 1,newBoundaryVerts[0], newBoundaryVerts[boundaryVerts.Count - 1], 
                secondExtrusionOffset + boundaryVerts.Count - 1, secondExtrusionOffset,  newBoundaryVerts[0],
        });


        int numSeeds = (int) ((Mathf.PerlinNoise1D(waterFallsSeed * 0.646385f) % 1.0f) * 3 + 2);
        int[] riverVerts = new int[numSeeds];

        int vert = (int)((Mathf.PerlinNoise1D(waterFallsSeed * 0.7432f) % 1.0f) * newBoundaryVerts.Count);

        riverVerts[0] = newBoundaryVerts[vert];

        for(int i = 1; i < numSeeds; i++) {
            riverVerts[i] = newBoundaryVerts[ ((vert * i) + (newBoundaryVerts.Count / numSeeds)) % newBoundaryVerts.Count];
        }

        for(int i = 1; i < numSeeds; i++) {
            Debug.DrawLine(verts[riverVerts[i-1]], verts[riverVerts[i]]);
        }
        Debug.DrawLine(verts[riverVerts[riverVerts.Length - 1]], verts[riverVerts[0]]);

        Vector3 avg = new();
        foreach (int item in riverVerts)
        {
            avg += verts[item];            
        }

        avg /= numSeeds;
        avg += new Vector3(Mathf.PerlinNoise1D(waterFallsSeed * 0.6382f), 0, Mathf.PerlinNoise1D(waterFallsSeed * 0.7685f)).normalized * 0.1f * radius;

        river_center = avg;

        for(int i = 0; i < numSeeds; i++) {
            Debug.DrawLine(verts[riverVerts[i]], avg, Color.red);
        }

        List<Vector3> riverCurves = new();

        for(int i = 0; i < numSeeds; i++) {
            int start = riverVerts[i];
            int end = riverVerts[(i + 1) % numSeeds];

            for(float t = 0; t < 1.0f; t += 0.1f) {
                Debug.DrawLine(EvaluateThreePointCurve(verts[start], verts[end], avg, t), 
                            EvaluateThreePointCurve(verts[start], verts[end], avg, t + 0.1f),
                            Color.green);
                
            }

            for(float t = 0; t <= 1.0f; t += 0.01f) {
                riverCurves.Add(EvaluateThreePointCurve(verts[start], verts[end], avg, t));
            }
        }


        for (int i = 0; i < terrainVerts.Count; i++) {
            float minDist = float.MaxValue;
            for(int j = 0; j < riverCurves.Count; j++) {
                float dist = (verts[terrainVerts[i]] - riverCurves[j]).magnitude;
                
                minDist = Mathf.Min(minDist, dist);
            }
            float distToCenter = (verts[terrainVerts[i]] - avg).magnitude;
            if (minDist < riverWidth && distToCenter > centerRiverCircleRadius) {
                Vector3 riverTranslation;
                if (verts[terrainVerts[i]].magnitude < 0.8f * radius) {
                    riverTranslation = new Vector3(0, 1f * (1 - minDist / riverWidth), 0);
                } else {
                    Vector3 temp = new Vector3(0, 1f * (1 - minDist / riverWidth), 0);
                    float t = (verts[terrainVerts[i]].magnitude - 0.8f * radius) / radius * 5f;
                    // float  t = 1;
                    riverTranslation = (1-t) * temp;
                }
                verts[terrainVerts[i]] -= riverTranslation;

                water.transform.position = new(
                    water.transform.position.x, 
                    verts[terrainVerts[i]].y,
                    water.transform.position.z
                );
            } else {
                verts[terrainVerts[i]] += new Vector3(0f, 0.5f*Mathf.PerlinNoise(terrainNoiseScale*verts[terrainVerts[i]].x, terrainNoiseScale*verts[terrainVerts[i]].z), 0f);
            }
        }


        water.transform.position -=  new Vector3(0, 1f, 0);


        mesh.vertices = verts.ToArray();
        mesh.triangles = tris.ToArray();
        mesh.RecalculateNormals();

        i++;
    }



    private Vector3 river_center;



    private float EvaluateY(float xPercent, float zPercent) {
        float normal = Evaluate2DNormal(xPercent, zPercent);

        float y;

        if(noise) {
            y = Mathf.PerlinNoise(xPercent * noiseScale, zPercent * noiseScale) * normal;
        }  else {
            y = normal;
        }

        return y;

    }




    public float stdDev = 0.2f;


    // Evaluate the 2D normal distri bution at given (x, y) coordinates
    public float Evaluate2DNormal(float x, float y)
    {
        // Mean values (centered at 0.5, 0.5 for normalized values)
        float meanX = 0.5f;
        float meanY = 0.5f;

        // Standard deviations (example: 0.1 for both x and y directions)
        float sigmaX = stdDev;
        float sigmaY = stdDev;

        // Assuming no correlation (rho = 0)
        float rho = 0.0f;

        // Calculate the normalization constant (1 / 2πσxσy * sqrt(1 - ρ^2))
        float normConst = 1.0f / (2.0f * Mathf.PI * sigmaX * sigmaY * Mathf.Sqrt(1.0f - rho * rho));

        // Compute the exponent part
        float exponent = -0.5f * (Mathf.Pow((x - meanX) / sigmaX, 2) +
                                 Mathf.Pow((y - meanY) / sigmaY, 2));

        // Return the evaluated PDF value
        return normConst * Mathf.Exp(exponent) / 12f;
    }



    public Vector3 EvaluateThreePointCurve(Vector3 start, Vector3 end, Vector3 control, float t) {
        Vector3 a = Mathf.Pow(1 - t, 2) * start;
        Vector3 b = 2 * (1 - t) * t * control;
        Vector3 c = Mathf.Pow(t, 2) * end;

        return a + b + c;
    }






    public void RenderWaterMesh() {
        List<Vector3> verts = new();
        List<int> tris = new();

        Dictionary<Vector3, int> verts_to_index = new();
        
        for(int i = -xRadius + 1; i < xRadius; i++) {
            for(int j = -zRadius + 1; j < zRadius; j++) {


                Vector3[] newVerts = {
                    new(i - 0.5f, 0, j - 0.5f),
                    new(i + 0.5f, 0, j - 0.5f),
                    new(i + 0.5f, 0, j + 0.5f),
                    new(i - 0.5f, 0, j + 0.5f),
                };

                int[] indicies = {
                    0, 0, 0, 0
                };


                int ii = 0;
                foreach (Vector3 vertex in newVerts) {
                    int index;
                    if (verts_to_index.ContainsKey(vertex)) {
                        index = verts_to_index[vertex];
                        indicies[ii] = index;
                    } else {
                        index = verts.Count;
                        verts.Add(vertex);

                        verts_to_index.Add(vertex, index);
                        indicies[ii] = index;
                    }

                    ii++;
                }

                //create triangles
                tris.AddRange(new List<int> {
                    indicies[0],  indicies[3], indicies[1],
                    indicies[1], indicies[3],  indicies[2],
                });
            }
        }


        List<Vector2> uvs = new();

        for(int i = 0; i < verts.Count; i++) {
            float u = xRadius / verts[i].x;
            float v = zRadius / verts[i].z;

            // x or z value is this percentage to the edge
            float xRatioToEdge = Math.Abs(u);
            float zRatioToEdge = Math.Abs(v);

            uvs.Add(new Vector2(1.0f/u + 1f, 1.0f/v + 1f) * 0.5f);

            float ratioToEdge = Mathf.Min(xRatioToEdge, zRatioToEdge);

            float length = (verts[i] * ratioToEdge).magnitude;

            float boundaryRatio = radius / length;

            // Debug.DrawLine(new Vector3(0, 0, 0), verts[i] * boundaryRatio);

            verts[i] *= boundaryRatio;
        }

        water.mesh = new Mesh();
        water.mesh.vertices = verts.ToArray();
        water.mesh.triangles = tris.ToArray();
        water.mesh.uv = uvs.ToArray();

        water.gameObject.GetComponent<Renderer>().material.SetVector(Shader.PropertyToID("center"), river_center);
    }
}
