import asyncio
import httpx
import os


async def get_size(client, url):
    response = await client.head(url=url)
    return int(response.headers["Content-Length"])


async def download_range(client, url, start, end, output):
    async with client.stream(
        "GET", url, headers={"Range": f"bytes={start}-{end}"}
    ) as response:
        with open(output, "wb") as f:
            async for chunk in response.aiter_raw():
                f.write(chunk)


async def download(url, output, chunks=16):
    async with httpx.AsyncClient() as client:
        file_size = await get_size(client, url)

        # create chunks amount of ranges based on file_size
        ranges = []
        if file_size > 10 * 1024 * 1024:
            for i in range(chunks):
                start = i * file_size // chunks
                end = min((i + 1) * file_size // chunks, file_size) - 1
                ranges.append((start, end, f"{output}.{i}"))
        else:
            ranges.append((0, file_size, output))

        await asyncio.gather(*[download_range(client, url, *range) for range in ranges])

        if len(ranges) == 1:
            return

        with open(output, "wb") as o:
            for _, _, chunk_path in ranges:
                with open(chunk_path, "rb") as s:
                    o.write(s.read())
                os.remove(chunk_path)


async def download_json(url):
    async with httpx.AsyncClient() as client:
        response = await client.get(url)
        return response.json()


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


async def setup_async(db="afdb_foldcomp", download_chunks=16):
    for i in ["", ".index", ".dbtype", ".lookup", ".source"]:
        await download(
            f"https://foldcomp.steineggerlab.workers.dev/{db}{i}",
            f"{db}{i}",
            chunks=download_chunks,
        )
