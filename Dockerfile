FROM python:3.12-slim
WORKDIR /dnabc
COPY . .
RUN pip install --no-cache-dir .
CMD ["bash"]