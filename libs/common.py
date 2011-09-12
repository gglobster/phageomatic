from os import path, makedirs

def ensure_dir(dir_path):
    """Check that the directory exists; if not, create it."""

    abs_path = path.abspath(dir_path)
    if not path.exists(abs_path):
        try: makedirs(abs_path)
        except Exception as message:
            status = 1
            # TODO: make graceful fail or request input if interactive mode
        else:
            message = 'created path'
            status = 0
    else:
        message = 'path exists'
        status = 0
    report = {'message': message, 'status': status}
    return abs_path, report

