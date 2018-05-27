﻿using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class TestRplyLoader : MonoBehaviour {
    public string fileName;
    public bool isAbsolute;
    public string absFilePath;

    GameObject go, go2;

    public Shader shader;

    int[] indices;

    List<GameObject> gos = new List<GameObject>();

    float curTime, prevTime;

    public int limitCount;

    void createMesh(int startIdx, int verticeCnt, ref Vector3[] vertex, ref Color32[] color)
    {
        Mesh mesh = new Mesh();
        //mesh.vertices = new Vector3[verticeCnt];

        Vector3[] curV = new Vector3[verticeCnt];
        Array.Copy(vertex, startIdx * limitCount, curV, 0, verticeCnt);
        mesh.vertices = curV;

        Color32[] curC = new Color32[verticeCnt];
        Array.Copy(color, startIdx * limitCount, curC, 0, verticeCnt);
        mesh.colors32 = curC;

        if (indices.Length > verticeCnt)
        {
            int[] subindices = new int[verticeCnt];
            Array.Copy(indices, subindices, verticeCnt);
            mesh.SetIndices(subindices, MeshTopology.Points, 0);
        }
        else
            mesh.SetIndices(indices, MeshTopology.Points, 0);
        mesh.name = "mesh" + startIdx.ToString();

        GameObject go = new GameObject("go" + startIdx.ToString());
        go.transform.parent = transform;

        MeshFilter mf = go.AddComponent<MeshFilter>();
        mf.mesh = mesh;

        MeshRenderer mr = go.AddComponent<MeshRenderer>();
        mr.material = new Material(shader);
    }


    void loadPLYDownSample()
    {
        if (!isAbsolute)
            absFilePath = System.IO.Path.Combine(Application.streamingAssetsPath, fileName);
        //zhenyi
        IntPtr plyIntPtr = PlyLoaderDll.LoadPly(absFilePath);
        
        Mesh mesh = new Mesh();
        Vector3[] vertices = PlyLoaderDll.GetRVertices(plyIntPtr);
        Color32[] colors = PlyLoaderDll.GetRColors(plyIntPtr);
        //int[] indices = PlyLoaderDll.GetRIndexs(plyIntPtr);
        PlyLoaderDll.UnLoadPly(plyIntPtr);

        int meshCount = vertices.Length / limitCount + 1;
        for (int i = 0; i < meshCount; i++)
        {
            createMesh(i, Math.Min(limitCount, vertices.Length - i * limitCount), ref vertices, ref colors);
        }
        
        
    }

    // Use this for initialization
    void Start () {
        indices = new int[limitCount];
        for (int i = 0; i < limitCount; i++)
        {
            indices[i] = i;
        }

        prevTime = Time.realtimeSinceStartup;
        loadPLYDownSample();
        curTime = Time.realtimeSinceStartup;
        print("whole took " + (curTime - prevTime) + "s");
    }
	
	// Update is called once per frame
	void Update () {
		
	}
}
