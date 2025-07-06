FROM python:3.10.19-slim

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Set working directory
WORKDIR /app

# Install dependencies
COPY requirements.txt .
RUN pip install --upgrade pip && pip install -r requirements.txt

# Copy project files
COPY . .

# Default command for the container (can be overridden)
CMD ["python", "run.py"]

# For deployment
# Run with Gunicorn (adjust module path if needed)
# CMD ["gunicorn", "-w", "4", "-b", "0.0.0.0:5000", "run:app"]