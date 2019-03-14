extern crate mazth;

use self::mazth::i_bound::IBound;

/// acceleration interface for building and querying spatial data
pub trait ISpatialAccel<T> where T: Default + Clone {
    /// query for a list of objects intersecting with input
    fn query_intersect( & self, input: &IBound ) -> Result< Vec< T >, & 'static str >;
    fn query_intersect_single( & self, input: &IBound ) -> Result< Vec< T >, & 'static str >;
    /// build a acceleration structure with input bounds and object ids
    fn build_all( & mut self, &[ (T, &IBound) ] ) -> Result< (), & 'static str >;
}
