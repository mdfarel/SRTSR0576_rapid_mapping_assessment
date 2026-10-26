// Cloudflare Worker for waterfields.geoai.id: serves webgis/ as static assets (see wrangler.jsonc).
// Static assets answer a Range request with the whole file, but PMTiles reads each archive in
// small byte ranges, so Range requests under /tiles/ are sliced here. Deploy with `make deploy`.
export default {
  async fetch(request, env) {
    const url = new URL(request.url);
    const range = request.headers.get("Range");
    if (!url.pathname.startsWith("/tiles/") || !range) return env.ASSETS.fetch(request);

    const whole = await env.ASSETS.fetch(new Request(url, { method: "GET" }));
    if (whole.status !== 200) return whole;
    const buf = await whole.arrayBuffer();
    const size = buf.byteLength;
    const m = /^bytes=(\d*)-(\d*)$/.exec(range.trim());
    if (!m || (m[1] === "" && m[2] === "")) return new Response(buf, whole);

    let start;
    let end;
    if (m[1] === "") {  // "bytes=-500": the last 500 bytes
      start = Math.max(0, size - Number(m[2]));
      end = size - 1;
    } else {
      start = Number(m[1]);
      end = m[2] === "" ? size - 1 : Math.min(Number(m[2]), size - 1);
    }
    if (start >= size || start > end) {
      return new Response(null, { status: 416, headers: { "Content-Range": `bytes */${size}` } });
    }
    const headers = new Headers(whole.headers);
    headers.set("Content-Range", `bytes ${start}-${end}/${size}`);
    headers.set("Content-Length", String(end - start + 1));
    headers.set("Accept-Ranges", "bytes");
    headers.set("Cache-Control", "public, max-age=3600");
    return new Response(buf.slice(start, end + 1), { status: 206, headers });
  },
};
