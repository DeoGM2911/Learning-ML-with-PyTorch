import sys
from pathlib import Path

services_path = str(Path(__file__).resolve().parent / "services")
if services_path not in sys.path:
    sys.path.append(services_path)

from app import create_app

app = create_app()

# Comment out the following line if you want to run the app with Gunicorn
if __name__ == "__main__":
    app.run(debug=True, use_reloader=False)