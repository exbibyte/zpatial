pub extern crate mazth;
extern crate pretty_env_logger;

pub use self::mazth::bound::{ AxisAlignedBBox, Axis };
pub use self::mazth::i_bound::{ IBound, BoundType };

use std::collections::HashMap;
use std::boxed::Box;
use std::f64;

use interface::i_spatial_accel::ISpatialAccel;
use interface::i_stat_tree::IStatTree;

/// implementation of spatial acceleration using bounding volume hierarchy with surface area heuristic
pub struct Bvh <T> where T: Default + Clone {
    _root: NodeBvh<T>,
    _bins: u32,
}

///internal node structure for Bvh
pub struct NodeBvh < T > where T: Default + Clone {
    _bound: AxisAlignedBBox,    
    _left: BvhBranch<T>,
    _right: BvhBranch<T>,
    _obj: T, //leaf data
}

pub enum BvhBranch<T> where T: Default + Clone {
    CHILD(Box< NodeBvh<T> >),
    EMPTY,
}

impl<T> Default for NodeBvh< T > where T: Default + Clone {
    fn default() -> NodeBvh<T> {
        NodeBvh {
            _bound: AxisAlignedBBox {
                _bound_lower: [ f64::NEG_INFINITY; 3 ],
                _bound_upper: [ f64::INFINITY; 3 ],
            },
            _left: BvhBranch::EMPTY,
            _right: BvhBranch::EMPTY,
            _obj: T::default(),
        }
    }
}

impl<T> NodeBvh < T > where T: Default + Clone {
    pub fn init_branches( b: AxisAlignedBBox, l: BvhBranch<T>, r: BvhBranch<T> ) -> NodeBvh<T> {
        NodeBvh {
            _bound: b,
            _left: l,
            _right: r,
            _obj: T::default(),
        }
    }
    pub fn init_leaf( b: AxisAlignedBBox, o: T ) -> NodeBvh<T> {
        NodeBvh {
            _bound: b,
            _left: BvhBranch::EMPTY,
            _right: BvhBranch::EMPTY,
            _obj: o,
        }
    }
    pub fn build_node( & mut self, num_bins: u32, objs: &[ (T, &IBound) ] ) -> Result< (), & 'static str > {
        for i in objs {
            match i.1.get_type() {
                BoundType::AxisAlignBox => (),
                _ => { return Err( "unsupported bound type" ) },
            }
        }
        let b = objs.iter().cloned().map( |x| x.1 ).collect::<Vec< &IBound > >();
        let mut u : AxisAlignedBBox = Default::default();
        u.get_union( &b[..] );

        // println!( "bound union: {:?}", u );

        //check for leaf condition
        if objs.len() == 1 {
            self._bound = u;
            self._obj = objs[0].0.clone();
            return Ok( () )
        }

        //todo: execute alternative branch when number of objects are low
        
        //get the longest axis of the bounding box
        let ( ref axis, ref length ) = u.get_longest_axis();

        self._bound = u;

        // // println!( "axis length: {}", length );
        
        let _c = self._bound.get_centroid();
        let bin_len = length / num_bins as f64;

        let mut bin_offsets = vec![];
        for i in 0..num_bins {
            match *axis {
                Axis::X => { bin_offsets.push( i as f64 * bin_len + self._bound._bound_lower[0] ); },
                Axis::Y => { bin_offsets.push( i as f64 * bin_len + self._bound._bound_lower[1] ); },
                Axis::Z => { bin_offsets.push( i as f64 * bin_len + self._bound._bound_lower[2] ); },
            }
        }

        //bin objects using use surface area heuristics with their entry and exit points in the long axis of bounding box
        //determine median and bin objects with their centroids
        let mut bins_surf_area = vec![ 0i64; num_bins as usize ];
        let mut obj_bin = vec![];
        let epsilon = 0.001f64;

        let mut centroids = vec![];
        
        for (idx,i) in objs.iter().enumerate() {
            let c = i.1.get_centroid();
            let obj_bound_data = i.1.get_bound_data();
            let obj_entry = &obj_bound_data[0..3];
            let obj_exit = &obj_bound_data[3..6];
            
            // println!("centroid: {:?}", c );
            let bin_id = match *axis {
                Axis::X => {
                    let entry_bin = ( ( obj_entry[0] - self._bound._bound_lower[0] ) / length * (num_bins as f64) * ( 1.0 - epsilon ) ) as usize;
                    let exit_bin = ( ( obj_exit[0] - self._bound._bound_lower[0] ) / length * (num_bins as f64) * ( 1.0 - epsilon ) ) as usize;
                    bins_surf_area[ entry_bin ] += 1;
                    bins_surf_area[ exit_bin ] -= 1;
                        
                    ( ( c[0] - self._bound._bound_lower[0] ) / length * (num_bins as f64) * ( 1.0 - epsilon ) ) as usize
                },
                Axis::Y => {
                    let entry_bin = ( ( obj_entry[1] - self._bound._bound_lower[1] ) / length * (num_bins as f64) * ( 1.0 - epsilon ) ) as usize;
                    let exit_bin = ( ( obj_exit[1] - self._bound._bound_lower[1] ) / length * (num_bins as f64) * ( 1.0 - epsilon ) ) as usize;
                    bins_surf_area[ entry_bin ] += 1;
                    bins_surf_area[ exit_bin ] -= 1;

                    ( ( c[1] - self._bound._bound_lower[1] ) / length * (num_bins as f64) * ( 1.0 - epsilon ) ) as usize
                },
                Axis::Z => {
                    let entry_bin = ( ( obj_entry[2] - self._bound._bound_lower[2] ) / length * (num_bins as f64) * ( 1.0 - epsilon ) ) as usize;
                    let exit_bin = ( ( obj_exit[2] - self._bound._bound_lower[2] ) / length * (num_bins as f64) * ( 1.0 - epsilon ) ) as usize;
                    bins_surf_area[ entry_bin ] += 1;
                    bins_surf_area[ exit_bin ] -= 1;

                    ( ( c[2] - self._bound._bound_lower[2] ) / length * (num_bins as f64) * ( 1.0 - epsilon ) ) as usize
                },
            };

            match *axis {
                Axis::X => {
                    centroids.push( (c[0],idx) );
                },
                Axis::Y => {
                    centroids.push( (c[1],idx) );
                },
                Axis::Z => {
                    centroids.push( (c[2],idx) );
                },
            }
            
            obj_bin.push( bin_id );
        }

        // println!( "bins surf area: {:?}", bins_surf_area );
        let mut bin_left : Vec< ( T, &IBound ) > = vec![];
        let mut bin_right : Vec< ( T, &IBound ) > = vec![];

        {   
            let _ = bins_surf_area.iter_mut().fold( 0, | acc, x | {
                *x += acc;
                *x
            }
            );

            let mut sum_surf_area = 0;
            for i in bins_surf_area.iter() {
                sum_surf_area += *i;
            }

            let median_surf_area = sum_surf_area / 2;
            // println!( "bins surf area: {:?}, sum_surf_area: {}", bins_surf_area, sum_surf_area );
            // println!( "median surf area: {:?}", median_surf_area );

            let mut split_bin_idx = 0;
            let mut accum_surf_area = 0;
            for (idx, i) in bins_surf_area.iter().enumerate() {
                accum_surf_area += *i;
                if accum_surf_area > median_surf_area {
                    split_bin_idx = idx;
                    break;
                }
            }
            
            let mut idx_obj = 0;
            for i in obj_bin {
                let obj = objs[ idx_obj ].clone();
                if i <= split_bin_idx {
                    bin_left.push( obj );
                } else {
                    bin_right.push( obj );
                }
                idx_obj += 1;
            }
        }
        
        if (bin_left.len() >= 2 && bin_right.len() == 0)||
           (bin_left.len() == 0 && bin_right.len() >=2) {
                bin_left.clear();
                bin_right.clear();
                centroids.sort_by(|a,b| a.partial_cmp(b).unwrap());
                let half = centroids.len()/2;
                for idx_ord in 0..half {
                    let obj = objs[ centroids[idx_ord].1 ].clone();
                    bin_left.push( obj );
                }
                for idx_ord in half..centroids.len() {
                    let obj = objs[ centroids[idx_ord].1 ].clone();
                    bin_right.push( obj );
                }      
        }
                
        // debug!( "split bin left count: {}, bin right count: {}", bin_left.len(), bin_right.len() );
        
        if bin_left.len() > 0 {
            let mut l : NodeBvh<T> = Default::default();
            // println!("num left children: {}", bin_left.len() );
            l.build_node( num_bins, &bin_left[..] )?;
            self._left = BvhBranch::CHILD( Box::new( l ) );
        } else {
            self._left = BvhBranch::EMPTY;
        }

        if bin_right.len() > 0 {
            let mut r : NodeBvh<T> = Default::default();
            // println!("num right children: {}", bin_right.len() );
            r.build_node( num_bins, &bin_right[..] )?;
            self._right = BvhBranch::CHILD( Box::new( r ) );
        } else {
            self._right = BvhBranch::EMPTY;
        }

        Ok( () )
    }
    pub fn search< F >( n : & NodeBvh<T>, b: & IBound, mut f : F ) where F : FnMut( T ) -> bool
    {
        let mut q = vec![ n ];
        while q.len() > 0 {
            let l = q.pop().unwrap();
            if l._bound.intersect( b ) {
                let mut present_l = true;
                let mut present_r = true;
                match l._left {
                    BvhBranch::CHILD( ref o ) => {
                        let o_ref : &NodeBvh<T> = o;
                        q.push( o_ref );
                    },
                    _ => present_l = false,
                }
                match l._right {
                    BvhBranch::CHILD( ref o ) => {
                        let o_ref : &NodeBvh<T> = o;
                        q.push( o_ref );
                    },
                    _ => present_r = false,
                }
                if !present_l && !present_r {
                    if f( l._obj.clone() ) {
                        return
                    }
                }
            }
        }
    }
}

impl <T> ISpatialAccel<T> for Bvh<T> where T: Default + Clone {
    fn query_intersect( & self, input: &IBound ) -> Result< Vec< T >, & 'static str >
    {
        match input.get_type() {
            BoundType::AxisAlignBox => (),
            _ => { return Err( "unsupported bound type" ) },
        }
        let mut out = vec![];
        {
            let func_collect = | x | { out.push( x ); false };
            NodeBvh::search( &self._root, input, func_collect );
        }
        Ok( out )
    }
    fn query_intersect_single( & self, input: &IBound ) -> Result< Vec< T >, & 'static str >
    {
        match input.get_type() {
            BoundType::AxisAlignBox => (),
            _ => { return Err( "unsupported bound type" ) },
        }
        let mut out = vec![];
        {
            let func_collect = | x | { out.push( x ); true };
            NodeBvh::search( &self._root, input, func_collect );
        }
        Ok( out )
    }
    fn build_all( & mut self, objs: &[ (T, &IBound) ] ) -> Result< (), & 'static str >
    {
        //initiate top down construction
        if self._bins == 0 {
            return Err( "bvh bin count cannot be zero" )
        }
        self._root.build_node( self._bins, objs )
    }
}

impl <T> Bvh<T> where T: Default + Clone {
    pub fn init( bins: u32 ) -> Bvh<T> {
        assert!( bins != 0 );
        Bvh {
            _root: NodeBvh {
                _bound: AxisAlignedBBox {
                    _bound_lower: [ f64::NEG_INFINITY; 3 ],
                    _bound_upper: [ f64::INFINITY; 3 ],
                },
                _left: BvhBranch::EMPTY,
                _right: BvhBranch::EMPTY,
                _obj: Default::default(),
            },
            _bins: bins,
        }
    }
}

impl<T> IStatTree for Bvh<T> where T: Default + Clone {
    fn sum_subtree_child_count( & self ) -> Option< u64 > {
        unimplemented!();
    }
    fn node_degree( & self ) -> Option< u64 > {
        unimplemented!();
    }
    fn max_depth( & self ) -> Option< u64 > {
        unimplemented!();
    }
    fn mean_depth( & self ) -> Option< u64 > {
        unimplemented!();
    }
    fn min_depth( & self ) -> Option< u64 > {
        unimplemented!();
    }
    fn balance_ratio( & self ) -> Option< f64 > {
        unimplemented!();
    }
}

