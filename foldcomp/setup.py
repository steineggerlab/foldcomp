import asyncio
import httpx
import os


async def get_size(client, url):
    response = await client.head(url=url)
    return int(response.headers["Content-Length"])


async def download_range(client, url, start, end, output, mode):
    curr_start = start
    with open(output, mode) as f:
        while True:
            try:
                async with client.stream(
                    "GET", url, headers={"Range": f"bytes={curr_start}-{end}"}
                ) as response:
                    async for chunk in response.aiter_raw():
                        f.write(chunk)
            except:
                # sometimes the connection is closed by the server
                # so we retry the download
                f.flush()
                os.fsync(f.fileno())
                # set start to current written position
                curr_start = start + f.tell()
                if curr_start < end:
                    continue
            break


async def download(url, output, chunks=16):
    async with httpx.AsyncClient() as client:
        file_size = await get_size(client, url)

        # check that file is not already downloaded
        if os.path.exists(output) and os.path.getsize(output) == file_size:
            return

        # create chunks amount of ranges based on file_size
        ranges = []
        if file_size > 10 * 1024 * 1024:
            for i in range(chunks):
                start = i * file_size // chunks
                end = min((i + 1) * file_size // chunks, file_size) - 1
                output_chunk = f"{output}.{i}"
                # adjust start if output_chunk already partially exists
                mode = "wb"
                if os.path.exists(output_chunk):
                    start += os.path.getsize(output_chunk)
                    mode = "ab"
                if start < end:
                    ranges.append((start, end, output_chunk, mode))
        else:
            ranges.append((0, file_size, output, "wb"))

        await asyncio.gather(*[download_range(client, url, *range) for range in ranges])

        if len(ranges) == 1:
            return

        with open(output, "wb") as o:
            BUFFER_SIZE = 10 * 1024 * 1024
            for _, _, chunk_path, _ in ranges:
                with open(chunk_path, "rb") as s:
                    while True:
                        buffer = s.read(BUFFER_SIZE)
                        if not buffer:
                            break
                        o.write(buffer)
                os.remove(chunk_path)


async def download_json(url):
    async with httpx.AsyncClient() as client:
        response = await client.get(url)
        return response.json()


# def list():
#     return asyncio.run(download_json("https://foldseek.steineggerlab.workers.dev/manifest.json"))


def setup(db="afdb_swissprot", download_chunks=16):
    for i in ["", ".index", ".dbtype", ".lookup", ".source"]:
        asyncio.run(
            download(
                f"https://foldcomp.steineggerlab.workers.dev/{db}{i}",
                f"{db}{i}",
                chunks=download_chunks,
            )
        )


async def setup_async(db="afdb_swissprot", download_chunks=16):
    for i in ["", ".index", ".dbtype", ".lookup", ".source"]:
        await download(
            f"https://foldcomp.steineggerlab.workers.dev/{db}{i}",
            f"{db}{i}",
            chunks=download_chunks,
        )
