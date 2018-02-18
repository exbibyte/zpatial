use interface::i_bound::IBound;

pub enum ClusterSize {
    Auto,
    Manual( usize ),
}

pub enum ClusterMetric {
    L1,
    L2,
}

pub enum ClusterMethod {
    Density,
    Average,
    Gaussian,
}

pub enum ClusterOrder {
    InOrder,
    PreSample,
}

/// acceleration interface for building and querying spatial data
pub trait ICluster {
    fn query_cluster( & self, input: u64 ) -> Result< u64, & 'static str >;
    fn build_all( & mut self, ClusterSize, ClusterMethod, ClusterMetric, ClusterOrder, input: &[ (u64, &IBound) ] ) -> Result< (), & 'static str >;
}
