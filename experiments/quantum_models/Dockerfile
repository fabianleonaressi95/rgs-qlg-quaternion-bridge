FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /app
COPY node.py /app/node.py

# The benchmark opens ordinary TCP sockets and needs no extra Linux capabilities.
USER 10001:10001
CMD ["python", "/app/node.py"]
