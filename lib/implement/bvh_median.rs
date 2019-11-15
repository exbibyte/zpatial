extern crate mazth;
extern crate pretty_env_logger;

use self::mazth::bound::{Axis, AxisAlignedBBox};
use self::mazth::i_bound::{BoundType, IBound};

use std::boxed::Box;
use std::f64;

use interface::i_spatial_accel::ISpatialAccel;
use interface::i_stat_tree::IStatTree;

/// implementation of spatial acceleration using bounding volume hierarchy with surface area heuristic
pub struct Bvh<T>
where
    T: Default + Clone,
{
    _root: NodeBvh<T>,
    _bins: u32,
}

///internal node structure for Bvh
pub struct NodeBvh<T>
where
    T: Default + Clone,
{
    _bound: AxisAlignedBBox,
    _left: BvhBranch<T>,
    _right: BvhBranch<T>,
    _obj: T, //leaf data
}

pub enum BvhBranch<T>
where
    T: Default + Clone,
{
    CHILD(Box<NodeBvh<T>>),
    EMPTY,
}

impl<T> Default for NodeBvh<T>
where
    T: Default + Clone,
{
    fn default() -> NodeBvh<T> {
        NodeBvh {
            _bound: AxisAlignedBBox {
                _bound_lower: [f64::NEG_INFINITY; 3],
                _bound_upper: [f64::INFINITY; 3],
            },
            _left: BvhBranch::EMPTY,
            _right: BvhBranch::EMPTY,
            _obj: T::default(),
        }
    }
}

impl<T> NodeBvh<T>
where
    T: Default + Clone,
{
    pub fn init_branches(b: AxisAlignedBBox, l: BvhBranch<T>, r: BvhBranch<T>) -> NodeBvh<T> {
        NodeBvh {
            _bound: b,
            _left: l,
            _right: r,
            _obj: T::default(),
        }
    }
    pub fn init_leaf(b: AxisAlignedBBox, o: T) -> NodeBvh<T> {
        NodeBvh {
            _bound: b,
            _left: BvhBranch::EMPTY,
            _right: BvhBranch::EMPTY,
            _obj: o,
        }
    }
    pub fn build_node(
        &mut self,
        num_bins: u32,
        objs: &[(T, &dyn IBound)],
    ) -> Result<(), &'static str> {
        for i in objs {
            match i.1.get_type() {
                BoundType::AxisAlignBox => (),
                _ => return Err("unsupported bound type"),
            }
        }
        let b = objs
            .iter()
            .cloned()
            .map(|x| x.1)
            .collect::<Vec<&dyn IBound>>();

        let mut u: AxisAlignedBBox = Default::default();
        u.get_union(&b[..]);

        //check for leaf condition
        if objs.len() == 1 {
            self._bound = u;
            self._obj = objs[0].0.clone();
            return Ok(());
        } else if objs.len() == 0 {
            return Ok(());
        }

        //get the longest axis of the bounding box
        let (ref axis, ref _length) = u.get_longest_axis();

        self._bound = u;

        let mut centroids = vec![];

        for (idx, i) in objs.iter().enumerate() {
            let c = i.1.get_centroid();
            match *axis {
                Axis::X => {
                    centroids.push((c[0], idx));
                }
                Axis::Y => {
                    centroids.push((c[1], idx));
                }
                Axis::Z => {
                    centroids.push((c[2], idx));
                }
            }
        }

        // println!( "bins surf area: {:?}", bins_surf_area );
        let mut bin_left: Vec<(T, &dyn IBound)> = vec![];
        let mut bin_right: Vec<(T, &dyn IBound)> = vec![];

        centroids.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        let half = centroids.len() / 2;
        assert!(centroids.len() > 1);
        for idx_ord in 0..half {
            let obj = objs[centroids[idx_ord].1].clone();
            bin_left.push(obj);
        }
        for idx_ord in half..centroids.len() {
            let obj = objs[centroids[idx_ord].1].clone();
            bin_right.push(obj);
        }

        // debug!( "split bin left count: {}, bin right count: {}", bin_left.len(), bin_right.len() );

        if bin_left.len() > 0 {
            let mut l: NodeBvh<T> = Default::default();
            // println!("num left children: {}", bin_left.len() );
            l.build_node(num_bins, &bin_left[..])?;
            self._left = BvhBranch::CHILD(Box::new(l));
        } else {
            self._left = BvhBranch::EMPTY;
        }

        if bin_right.len() > 0 {
            let mut r: NodeBvh<T> = Default::default();
            // println!("num right children: {}", bin_right.len() );
            r.build_node(num_bins, &bin_right[..])?;
            self._right = BvhBranch::CHILD(Box::new(r));
        } else {
            self._right = BvhBranch::EMPTY;
        }

        Ok(())
    }
    pub fn search<F>(n: &NodeBvh<T>, b: &dyn IBound, mut f: F)
    where
        F: FnMut(T) -> bool,
    {
        let mut q = vec![n];
        while q.len() > 0 {
            let l = q.pop().unwrap();
            if l._bound.intersect(b) {
                let mut present_l = true;
                let mut present_r = true;
                match l._left {
                    BvhBranch::CHILD(ref o) => {
                        let o_ref: &NodeBvh<T> = o;
                        q.push(o_ref);
                    }
                    _ => present_l = false,
                }
                match l._right {
                    BvhBranch::CHILD(ref o) => {
                        let o_ref: &NodeBvh<T> = o;
                        q.push(o_ref);
                    }
                    _ => present_r = false,
                }
                if !present_l && !present_r {
                    if f(l._obj.clone()) {
                        return;
                    }
                }
            }
        }
    }
}

impl<T> ISpatialAccel<T> for Bvh<T>
where
    T: Default + Clone,
{
    fn query_intersect(&self, input: &dyn IBound) -> Result<Vec<T>, &'static str> {
        match input.get_type() {
            BoundType::AxisAlignBox => (),
            _ => return Err("unsupported bound type"),
        }
        let mut out = vec![];
        {
            let func_collect = |x| {
                out.push(x);
                false
            };
            NodeBvh::search(&self._root, input, func_collect);
        }
        Ok(out)
    }
    fn query_intersect_single(&self, input: &dyn IBound) -> Result<Vec<T>, &'static str> {
        match input.get_type() {
            BoundType::AxisAlignBox => (),
            _ => return Err("unsupported bound type"),
        }
        let mut out = vec![];
        {
            let func_collect = |x| {
                out.push(x);
                true
            };
            NodeBvh::search(&self._root, input, func_collect);
        }
        Ok(out)
    }
    fn build_all(&mut self, objs: &[(T, &dyn IBound)]) -> Result<(), &'static str> {
        //initiate top down construction
        if self._bins == 0 {
            return Err("bvh bin count cannot be zero");
        }
        self._root.build_node(self._bins, objs)
    }
}

impl<T> Bvh<T>
where
    T: Default + Clone,
{
    pub fn init(bins: u32) -> Bvh<T> {
        assert!(bins != 0);
        Bvh {
            _root: NodeBvh {
                _bound: AxisAlignedBBox {
                    _bound_lower: [f64::NEG_INFINITY; 3],
                    _bound_upper: [f64::INFINITY; 3],
                },
                _left: BvhBranch::EMPTY,
                _right: BvhBranch::EMPTY,
                _obj: Default::default(),
            },
            _bins: bins,
        }
    }
}

impl<T> IStatTree for Bvh<T>
where
    T: Default + Clone,
{
    fn sum_subtree_child_count(&self) -> Option<u64> {
        unimplemented!();
    }
    fn node_degree(&self) -> Option<u64> {
        unimplemented!();
    }
    fn max_depth(&self) -> Option<u64> {
        unimplemented!();
    }
    fn mean_depth(&self) -> Option<u64> {
        unimplemented!();
    }
    fn min_depth(&self) -> Option<u64> {
        unimplemented!();
    }
    fn balance_ratio(&self) -> Option<f64> {
        unimplemented!();
    }
}
