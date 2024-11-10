
using System;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UIElements;

public class VoxelCircle : MonoBehaviour
{
    private MeshFilter meshFilter;
    private Mesh mesh;

    public int xRadius = 40, zRadius = 40;


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

        for(int i = -xRadius + 1; i < xRadius; i++) {
            for(int j = -zRadius + 1; j < zRadius; j++) {

                float angle = Mathf.Atan(((float) i) / j);

                float radX = xRadius * Mathf.Sin(angle);
                float radZ = zRadius * Mathf.Cos(angle);

                if (Mathf.Abs(i) > Mathf.Abs(radX) && Mathf.Abs(j) > Mathf.Abs(radZ)) {
                    continue;
                }

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
                    if (verts_to_index.ContainsKey(vertex)) {
                        int index = verts_to_index[vertex];

                        indicies[ii] = index;
                    } else {
                        int index = verts.Count;
                        verts.Add(vertex);

                        verts_to_index.Add(vertex, index);
                        indicies[ii] = index;
                    }

                    ii++;
                }

                //create triangles
                tris.AddRange(new List<int> {
                    indicies[0], indicies[3], indicies[1],
                    indicies[1], indicies[3], indicies[2],
                });
            }
        }

        mesh.vertices = verts.ToArray();
        mesh.triangles = tris.ToArray();
        mesh.RecalculateNormals();
    }
}
