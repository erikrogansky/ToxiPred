import hashlib, json

def make_etag(payload: dict, last_modified: float) -> str:
  body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
  h = hashlib.sha256((body + str(last_modified)).encode("utf-8")).hexdigest()
  return f'W/"{h[:32]}"'
