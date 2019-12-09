FROM python:3

WORKDIR /app

COPY Docker/requirements.txt /app

COPY bin/ /app
COPY Diagnostic/ /app

RUN pip3 install --no-cache-dir -r requirements.txt
RUN apt-get update && apt-get install --no-install-recommends -qy gosu python3-tk

ADD . /app

ENTRYPOINT ["bash", "Docker/entrypoint.sh"]


