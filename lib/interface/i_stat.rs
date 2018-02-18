pub trait IStatRate {
    fn reset( & self ) -> Result< (), & 'static str >;
    fn avg( & self ) -> Option< ( f64/*period*/, f64/*rate*/ ) >;
    fn max( & self ) -> Option< ( f64, f64 ) >;
    fn min( & self ) -> Option< ( f64, f64 ) >;
    fn set_window_period( & self, f64 ) -> Result< (), & 'static str >;
}
