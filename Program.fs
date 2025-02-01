open System
open System.Numerics
open System.IO
open System.Diagnostics
open System.Threading.Tasks
open System.Threading

type Ray =
    struct
        val origin: Vector3
        val direction: Vector3

        new(origin, direction) =
            { origin = origin
              direction = direction }
    end

[<Struct>]
type Material =
    | Diffuse = 0
    | Specular = 1
    | Refractive = 2

type Sphere =
    struct
        val center: Vector3
        val radius: float32
        val emission: Vector3
        val material: Material
        val color: Vector3

        new(radius, center, emission, color, material) =
            { radius = radius
              center = center
              emission = emission
              color = color
              material = material }

        member inline this.intersect(ray: inref<Ray>) : float32 =
            let eps = 1e-3f
            let f = ray.origin - this.center
            let a = ray.direction.LengthSquared()
            let b = - Vector3.Dot(f, ray.direction)
            let c = f.LengthSquared() - this.radius * this.radius
            let delta = this.radius * this.radius - (f + b / a * ray.direction).LengthSquared()
            if delta < 0.0f then
                0.0f
            else
                let q = b + float32 (sign b) * sqrt (a * delta)
                let t0, t1 = c / q, q / a
                if t0 > eps then t0
                else if t1 > eps then t1
                else 0.0f

    end

let spheres =
    [| Sphere(1e4f, Vector3(1e4f + 1f, 40.8f, 81.6f), Vector3.Zero, Vector3(0.75f, 0.25f, 0.25f), Material.Diffuse)
       Sphere(1e4f, Vector3(-1e4f + 99f, 40.8f, 81.6f), Vector3.Zero, Vector3(0.25f, 0.25f, 0.75f), Material.Diffuse)
       Sphere(1e4f, Vector3(50f, 40.8f, 1e4f), Vector3.Zero, Vector3(0.75f, 0.75f, 0.75f), Material.Diffuse)
       Sphere(1e4f, Vector3(50f, 40.8f, -1e4f + 181f), Vector3.Zero, Vector3.Zero, Material.Diffuse)
       Sphere(1e4f, Vector3(50f, 1e4f, 81.6f), Vector3.Zero, Vector3(0.75f, 0.75f, 0.75f), Material.Diffuse)
       Sphere(1e4f, Vector3(50f, -1e4f + 81.6f, 81.6f), Vector3.Zero, Vector3(0.75f, 0.75f, 0.75f), Material.Diffuse)
       Sphere(16.5f, Vector3(27f, 16.5f, 47f), Vector3.Zero, Vector3(0.999f, 0.999f, 0.999f), Material.Specular)
       Sphere(16.5f, Vector3(73f, 16.5f, 78f), Vector3.Zero, Vector3(0.999f, 0.999f, 0.999f), Material.Refractive)
       Sphere(600f, Vector3(50f, 681.6f - 0.27f, 81.6f), Vector3(12f, 12f, 12f), Vector3.Zero, Material.Diffuse) |]

let inline intersect (ray: inref<Ray>) (t: byref<float32>) (id: byref<int32>) : bool =
    let inf = Single.PositiveInfinity
    t <- inf
    for i = 0 to spheres.Length - 1 do
        let d = spheres[i].intersect &ray
        if d <> 0f && d < t then
            t <- d
            id <- i
    t <> inf

let radiance (ray: inref<Ray>) (rng: Random) : Vector3 =
    let maxDepth = 8
    let eps = 1e-2f
    let mutable t = 0f
    let mutable id = 0
    let mutable depth = 0
    let mutable ray = ray
    let mutable f = Vector3.One
    let mutable r = Vector3.Zero
    while depth <= maxDepth do
        depth <- depth + 1
        if not (intersect &ray &t &id) then
            depth <- maxDepth + 1
        else
            let obj = spheres[id]
            let x = ray.origin + ray.direction * t
            let n = Vector3.Normalize(x - obj.center)
            let nl = if Vector3.Dot(n, ray.direction) < 0.0f then n else -n
            r <- r + f * obj.emission
            f <- f * obj.color
            let p = max f.X (max f.Y f.Z)
            let rrDepth = 5
            if depth > rrDepth && rng.NextSingle() >= p then
                depth <- maxDepth + 1
            else
                f <- if depth > rrDepth then f / p else f
                match obj.material with
                | Material.Diffuse ->
                    let r1 = 2.0f * MathF.PI * rng.NextSingle()
                    let r2 = rng.NextSingle()
                    let r2s = sqrt r2
                    let w = nl
                    let u =
                        if abs w.X > 0.1f then
                            Vector3(0f, 1f, 0f)
                        else
                            Vector3(1f, 0f, 0f)
                        |> fun x -> Vector3.Cross(x, w)
                        |> Vector3.Normalize
                    let v = Vector3.Cross(w, u)
                    let d = u * (cos r1) * r2s + v * (sin r1) * r2s + w * sqrt (1.0f - r2)
                    ray <- Ray(x + eps * nl, d)
                | Material.Specular ->
                    let d = ray.direction - n * 2f * Vector3.Dot(n, ray.direction)
                    ray <- Ray(x + eps * nl, d)
                | _ ->
                    let rdir = ray.direction - n * 2f * Vector3.Dot(n, ray.direction)
                    let reflectRay = Ray(x + eps * nl, rdir)
                    let into = Vector3.Dot(n, nl) > 0.0f
                    let nc, nt = 1.0f, 1.5f
                    let nnt = if into then nc / nt else nt / nc
                    let ddn = Vector3.Dot(ray.direction, nl)
                    let cos2t = 1.0f - nnt * nnt * (1.0f - ddn * ddn)
                    if cos2t < 0f then
                        ray <- reflectRay
                    else
                        let tdir =
                            (ray.direction * nnt) - (n * (if into then 1.0f else -1.0f) * (ddn * nnt + sqrt cos2t))
                        let refractRay = Ray(x - eps * nl, tdir)
                        let a, b = nt - nc, nt + nc
                        let R0, c = a * a / (b * b), 1f - (if into then -ddn else Vector3.Dot(tdir, n))
                        let Re = R0 + (1.0f - R0) * c * c * c * c * c
                        let Tr = 1.0f - Re
                        let P = 0.25f + 0.5f * Re
                        if rng.NextSingle() < P then
                            ray <- reflectRay
                            f <- f * Re / P
                        else
                            ray <- refractRay
                            f <- f * Tr / (1.0f - P)
    r

let inline toInt (x: float32) =
    int (MathF.Pow(x, 1f / 2.2f) * 255.0f + 0.5f)

[<EntryPoint>]
let main (argv: string[]) =
    let w, h, samples = 1024, 768, if argv.Length = 1 then int argv[0] else 64
    let c = Array.zeroCreate<Vector3> (w * h)
    let inline renderRow (y: int) (rng: Random) =
        let cam =
            Ray(Vector3(50f, 52f, 295.6f), Vector3.Normalize(Vector3(0f, -0.042612f, -1f)))
        let cx = Vector3(float32 w * 0.5135f / float32 h, 0.0f, 0.0f)
        let cy = Vector3.Normalize(Vector3.Cross(cx, cam.direction)) * 0.5135f
        for x = 0 to w - 1 do
            for sy = 0 to 1 do
                for sx = 0 to 1 do
                    let mutable r = Vector3.Zero
                    for s = 0 to samples - 1 do
                        let r1 = 2.0f * rng.NextSingle()
                        let dx =
                            if r1 < 1.0f then
                                (sqrt r1) - 1.0f
                            else
                                1.0f - (sqrt (2.0f - r1))

                        let r2 = 2.0f * rng.NextSingle()

                        let dy =
                            if r2 < 1.0f then
                                (sqrt r2) - 1.0f
                            else
                                1.0f - (sqrt (2.0f - r2))

                        let d =
                            cam.direction
                            + cx * (((float32 sx + 0.5f + dx) / 2.0f + float32 x) / float32 w - 0.5f)
                            + cy * (((float32 sy + 0.5f + dy) / 2.0f + float32 y) / float32 h - 0.5f)
                            |> Vector3.Normalize

                        let ray = Ray(cam.origin + 140f * d, d)
                        r <- r + (radiance &ray rng) / float32 samples
                    let i = x + (h - 1 - y) * w
                    c[i] <- c[i] + 0.25f * Vector3.Clamp(r, Vector3.Zero, Vector3.One)
    let stopwatch = Stopwatch.StartNew()
    seq {
        for y = 0 to h - 1 do
            let random = new ThreadLocal<Random>(fun () -> Random(y))
            yield Task.Run(fun () -> renderRow y random.Value)
    }
    |> Seq.toArray
    |> Task.WhenAll
    |> Async.AwaitTask
    |> Async.RunSynchronously
    stopwatch.Stop()
    printfn $"Elapsed time: %A{stopwatch.Elapsed}"
    use f = File.CreateText "image.ppm"
    f.Write $"P3\n{w} {h}\n{255}\n"
    for i = 0 to c.Length - 1 do
        f.Write $"{toInt c[i].X} {toInt c[i].Y} {toInt c[i].Z} "
    0
