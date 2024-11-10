using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;

public class TerrainRenderer : MonoBehaviour
{

    private MeshFilter meshFilter;
    private Mesh mesh;

    public float height = 50f;

    public int radius;

    public float voronoiScale = 1f;
    public float voronoiSeed = 534.234f;

    public float max;
    public float min;



    // Start is called before the first frame update
    void Start()
    {
        meshFilter = GetComponent<MeshFilter>();
        mesh = new();
    }

    // Update is called once per frame
    void Update()
    {
        RenderMesh();

        meshFilter.mesh = mesh;
    }




    void RenderMesh() {
        List<Vector3> verts = new();
        List<int> tris = new();

        Dictionary<Vector3, int> verts_to_index = new();


        for(int i = -radius + 1; i < radius; i++) {
            for(int j = -radius + 1; j < radius; j++) {

                float xPercent = (i - (-radius + 1)) / (radius - (-radius + 1)); 
                float zPercent = (j - (-radius + 1)) / (radius - (-radius + 1)); 

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

                    // float value = 20*Terraformer.Evaluate2DNormal((voronoi * vertex.x / radius) / 2 + 0.5f, (voronoi * vertex.z / radius) / 2 + 0.5f) *StackedPerlin(vertex.x / radius, vertex.z /radius, voronoiSeed);
                    float value = 0f;
                    // verts[index] += new Vector3(0, 0.2f*Voronoi.sampleVoronoi(vertex.x, vertex.z, voronoiScale, voronoiSeed, (int)(voronoiScale / 4)), 0);

                    ii++;
                }

                //create triangles
                tris.AddRange(new List<int> {
                    indicies[0], indicies[3], indicies[1], 
                    indicies[1],  indicies[3], indicies[2],
                });
            }
        }

        for (int i = 0; i < verts.Count; i++) {
            float voronoi = Voronoi.sampleVoronoi(verts[i].x, verts[i].z, voronoiScale, voronoiSeed, (int)(voronoiScale / 4));
            float value = 
                Evaluate2DNormal(
                    (verts[i].x / radius )/ 2 + 0.5f,
                    (verts[i].z / radius ) / 2 + 0.5f
                ) *
                2 * (StackedPerlin(verts[i].x / radius + voronoi, verts[i].z /radius + voronoi, voronoiSeed) - 0.5f);

            min = Mathf.Min(min, value);
            max = Mathf.Min(min, value);

            verts[i] += new Vector3(0, height * value, 0);
        }
         
         

        mesh.vertices = verts.ToArray();
        mesh.triangles = tris.ToArray();
        mesh.RecalculateNormals();
    }



    float StackedPerlin(float x, float y, float seed) {
        float first = Mathf.PerlinNoise(x + 0.75f*seed, y + 3.74f*seed);
        float second = Mathf.PerlinNoise(2*x + 547.54f*seed, 2*y + .37f*seed);
        float third = Mathf.PerlinNoise(4*x + 56.23f*seed, 4*y + 73.23f*seed);
        float fourth = Mathf.PerlinNoise(8*x + 56.65f*seed, 8*y + 95.11f*seed);

        return (first + second/2 + third/4 + fourth/8) / (1f + 0.5f + 0.25f + 0.125f);
    }


    
    [Range(0, 1)]
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
    
}
