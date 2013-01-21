"""Command-line interface for maintenance of the cgp package."""

if __name__ == "__main__":
    
    import argparse
    import sys
    
    parser = argparse.ArgumentParser(description=__doc__.strip(), 
        prog="python -m cgp")
    parser.add_argument("--clear-urlcache", action="store_true", 
        help="Clear web service cache")
    parser.add_argument("--clear-cellml2py", action="store_true", 
        help="Clear cache of Python modules autogenerated from CellML")
    parser.add_argument("--clear", action="store_true", 
        help="Clear all caches")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()
    if args.clear:
        args.clear_urlcache = True
        args.clear_cellml2py = True
    if args.clear_urlcache:
        from cgp.physmod.cellmlmodel import mem
        mem.clear()
    if args.clear_cellml2py:
        import shutil
        from cgp.physmod.cellmlmodel import cgp_tempdir
        print "Clearing " + cgp_tempdir
        shutil.rmtree(cgp_tempdir)
