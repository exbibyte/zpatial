extern crate chrono;
extern crate mazth;
extern crate rand;

use self::chrono::Local;
use self::rand::Rng;

use self::mazth::{
    bound::AxisAlignedBBox, bound_sphere::BoundSphere, i_bound::IBound, i_shape::IShape,
    i_shape::ShapeType, point::Point3, triprism::TriPrism,
};
use implement::bvh_median::Bvh;
use interface::i_spatial_accel::ISpatialAccel;
#[cfg(test)]
use std::f64;

#[test]
fn test_bvh_median_supported_bounds() {
    let mut a = Bvh::init(10); //bin count of 10
    let aabb = AxisAlignedBBox::init(ShapeType::Sphere, &[0f64, 0f64, 0f64, 5f64]);
    let objs = [(0u64, &aabb as &dyn IBound)];
    match a.build_all(&objs[..]) {
        Ok(()) => (),
        _ => {
            panic!("unexpected result for supported bound type");
        }
    }
}

#[test]
fn test_bvh_median_unsupported_bounds() {
    let mut a = Bvh::init(10);
    let aabb = BoundSphere::init(ShapeType::Sphere, &[0f64, 0f64, 0f64, 5f64]);
    let objs = [(0u64, &aabb as &dyn IBound)];
    match a.build_all(&objs[..]) {
        Err(_) => (),
        _ => {
            panic!("unexpected result for unsupported bound type");
        }
    }
}

#[test]
fn test_bvh_median_construction_and_query() {
    let mut a = Bvh::init(30);
    let mut bounds = vec![];
    let mut bound_refs = vec![];
    for i in 0..20 {
        let aabb = AxisAlignedBBox::init(
            ShapeType::Sphere,
            &[f64::from(i), f64::from(i), f64::from(i), 5f64],
        );
        bounds.push(aabb);
    }
    for i in 0..20 {
        bound_refs.push((i as u64, &bounds[i] as &dyn IBound));
    }

    match a.build_all(&bound_refs[..]) {
        Ok(()) => (),
        _ => {
            panic!("unexpected result for supported bound type");
        }
    }

    //query within
    {
        let i = 0f64;
        let query = AxisAlignedBBox::init(ShapeType::Point, &[i, i, i]);
        match a.query_intersect(&query) {
            Ok(o) => {
                assert!(o.len() == 6, "bvh query_intersect return length unexpected");
                for j in o {
                    assert!(j <= 5, "bvh query_intersect return index unexpected");
                }
            }
            _ => (),
        }
    }
    {
        let i = 19f64;
        let query = AxisAlignedBBox::init(ShapeType::Point, &[i, i, i]);
        match a.query_intersect(&query) {
            Ok(o) => {
                assert!(o.len() == 6, "bvh query_intersect return length unexpected");
                for j in o {
                    assert!(j >= 14, "bvh query_intersect return index unexpected");
                }
            }
            _ => (),
        }
    }
    {
        let i = -5f64;
        let query = AxisAlignedBBox::init(ShapeType::Point, &[i, i, i]);
        match a.query_intersect(&query) {
            Ok(o) => {
                assert!(o.len() == 1, "bvh query_intersect return length unexpected");
                for j in o {
                    assert!(j == 0, "bvh query_intersect return index unexpected");
                }
            }
            _ => (),
        }
    }

    //query not present
    {
        let i = -5.1f64;
        let query = AxisAlignedBBox::init(ShapeType::Point, &[i, i, i]);
        match a.query_intersect(&query) {
            Ok(o) => {
                assert!(o.len() == 0, "bvh query_intersect return length unexpected");
            }
            _ => (),
        }
    }
    {
        let i = 26f64;
        let query = AxisAlignedBBox::init(ShapeType::Point, &[i, i, i]);
        match a.query_intersect(&query) {
            Ok(o) => {
                assert!(o.len() == 0, "bvh query_intersect return length unexpected");
            }
            _ => (),
        }
    }
}

#[test]
fn test_bvh_median_stress() {
    // env::set_var("LOG_SETTING", "info" );

    // pretty_env_logger::init_custom_env( "LOG_SETTING" );

    let mut a = Bvh::init(10); //bin count of 10

    let mut rng = rand::thread_rng();

    let v = (0..100_000)
        .map(|x| {
            let rx = rng.gen_range(0., 1.);
            let ry = rng.gen_range(0., 1.);
            let rz = rng.gen_range(0., 1.);
            (
                x,
                AxisAlignedBBox::init(ShapeType::Sphere, &[rx, ry, rz, 0.000001f64]),
            )
        })
        .collect::<Vec<_>>();

    let objs = v
        .iter()
        .map(|x| (x.0, &x.1 as &dyn IBound))
        .collect::<Vec<_>>();

    let t0 = Local::now();

    info!("build tree");

    match a.build_all(&objs[..]) {
        Ok(()) => (),
        _ => {
            panic!("unexpected result for supported bound type");
        }
    }

    let t1 = Local::now();
    let t_delta = t1.signed_duration_since(t0).num_microseconds().unwrap() as f64;

    info!("time to build tree: {} microseconds", t_delta);

    let mut query_ret_count = 0;
    let mut query_time = 0.;

    for (idx, vv) in v.iter().enumerate() {
        let c = vv.1.get_centroid();
        let query = AxisAlignedBBox::init(ShapeType::Point, &c[..]);

        let t2 = Local::now();

        match a.query_intersect_single(&query) {
            Ok(o) => {
                // assert!( o.len() > 0, "bvh query_intersect return length unexpected" );

                let t3 = Local::now();

                let t_delta = t3.signed_duration_since(t2).num_microseconds().unwrap() as f64;
                query_time += t_delta;

                let found = o.iter().any(|x| *x as usize == idx);
                assert!(found, "query item not found");
                query_ret_count += o.len();
            }
            _ => {
                panic!("query unexpected result");
            }
        }
    }
    info!(
        "avg query return size: {}",
        query_ret_count as f32 / v.len() as f32
    );
    info!("avg query time: {}", query_time as f32 / v.len() as f32);
}

#[test]
fn test_bvh_median_tri_prism_point_intersection() {
    let mut a = Bvh::init(30);

    let triprism_base_verts = [0., 0., 0., 1., 0., 0., 1., 1., 0.];

    let tp = TriPrism::init(&triprism_base_verts, 1.);

    let mut bounds = vec![];

    let aabb = tp._bound.clone();
    bounds.push(aabb);

    let objs = &[(tp.clone(), &bounds[0] as &dyn IBound)];

    match a.build_all(&objs[..]) {
        Ok(()) => (),
        _ => {
            panic!("unexpected result for supported bound type");
        }
    }

    //query present
    {
        let b = Point3::init(&[0.5, 0.5, 0.5]);
        let query = b._bound.clone();

        match a.query_intersect(&query) {
            Ok(o) => {
                assert!(o.len() > 0, "bvh query_intersect return length unexpected");
                let intersected = o.iter().any(|x| x.get_intersect(&b).0);
                assert_eq!(intersected, true);
            }
            _ => (),
        }
    }
    //query not present
    {
        let b = Point3::init(&[0.5, 0.51, 0.5]);
        let query = b._bound.clone();

        match a.query_intersect(&query) {
            Ok(o) => {
                let intersected = o.iter().any(|x| x.get_intersect(&b).0);
                assert_eq!(intersected, false);
            }
            _ => (),
        }
    }
}
