import asyncio
import aiohttp
import os


async def get_size(session, url):
    async with session.head(url=url) as response:
        return int(response.headers["Content-Length"])


async def download_range(session, url, range):
    (start, end, output) = range
    async with session.get(
        url=url, headers={"Range": f"bytes={start}-{end}"}
    ) as response:
        with open(output, "wb") as f:
            # write StreamReader response.content to file f
            while True:
                chunk = await response.content.read(1024)
                if not chunk:
                    break
                f.write(chunk)


async def download(url, output, chunks=16):
    async with aiohttp.ClientSession(auto_decompress=False) as session:
        file_size = await get_size(session, url)

        # create N ranges based on file_size
        ranges = []
        for i in range(chunks):
            start = i * file_size // chunks
            end = min((i + 1) * file_size // chunks, file_size) - 1
            ranges.append((start, end, f"{output}.{i}"))

        await asyncio.gather(*[download_range(session, url, range) for range in ranges])

        with open(output, "wb") as o:
            for _, _, _, chunk_path in ranges:
                with open(chunk_path, "rb") as s:
                    o.write(s.read())
                os.remove(chunk_path)


async def download_json(url):
    async with aiohttp.ClientSession() as session:
        async with session.get(url) as response:
            return await response.json()


# def list():
#     return asyncio.run(download_json("https://foldseek.steineggerlab.workers.dev/manifest.json"))


def setup(db="afdb_foldcomp", download_chunks=16):
    for i in ["", ".index", ".dbtype", ".lookup", ".source"]:
        asyncio.run(
            download(
                f"https://foldcomp.steineggerlab.workers.dev/{db}{i}",
                f"{db}{i}",
                chunks=download_chunks,
            )
        )
