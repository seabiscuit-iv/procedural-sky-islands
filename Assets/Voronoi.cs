using System;
using Unity.VisualScripting;

public class Voronoi
{
    // Function to generate a deterministic "random" value for a given seed point
    private static float GenerateSeedValue(float seed, int index)
    {
        // A simple deterministic pseudo-random generator using seed and index
        return (float)(Math.Sin(seed + index) * 10000) % 1.0f;
    }

    // Function to generate a "random" value associated with the closest seed point
    public static float sampleVoronoi(float x, float y, float scale, float seed, int numSeeds)
    {
        // Define the number of seeds
        var seeds = new (float X, float Y)[numSeeds];

        // Generate seeds based on the seed value (for example, use the seed value to create 4 seed points)
        for (int i = 0; i < numSeeds; i++)
        {
            seeds[i].X = GenerateSeedValue(seed * 5.235f, i) * scale; // Scale seed value to control the seed locations
            seeds[i].Y = GenerateSeedValue(seed * 36.97f, i + 1) * scale; // Slight offset for variation'
        }

        // Now, compute the closest seed point to (x, y)
        int closestSeedIndex = -1;
        float minDistanceSquared = float.MaxValue;

        for (int i = 0; i < seeds.Length; i++)
        {
            float dx = x - seeds[i].X;
            float dy = y - seeds[i].Y;
            float distSquared = dx * dx + dy * dy; // squared distance to avoid sqrt calculation

            if (distSquared < minDistanceSquared)
            {
                minDistanceSquared = distSquared;
                closestSeedIndex = i;
            }
        }

        // Generate a "random" value for the closest seed
        // This is a deterministic value associated with that seed
        return GenerateSeedValue(seed, closestSeedIndex);
    }

}
