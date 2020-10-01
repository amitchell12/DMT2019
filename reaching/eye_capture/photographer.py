import codecs
import logging
import socket
import sys
import threading
import time
import traceback

import cv2

capture = None


class VCExhaustBuffer:
    """
    Exhaust the buffer every time, take the frame that takes longer than
    a time limit
    """

    def __init__(self, cam_id):
        self.cap = cv2.VideoCapture(cam_id)
        self.max_exhaust = 10  # prevent infinite loop
        self.time_threshold_ms = 5

    def set(self, *args, **kwargs):
        return self.cap.set(*args, **kwargs)

    def isOpened(self):
        return self.cap.isOpened()

    def read(self):
        # take last in buffer
        last = None
        for idx in range(self.max_exhaust):
            t0 = time.time()
            ok, im = self.cap.read()
            if not ok:
                return False, None
            else:
                delay_ms = 1000*(time.time() - t0)
                if delay_ms > self.time_threshold_ms:
                    return ok, im
            last = im
        return bool(last), last


class VCTimeThread:
    """
    Take the last frame that didn't exhaust the buffer, or if none, wait for
    a frame
    """
    def __init__(self, cam_id):
        self.cap = cv2.VideoCapture(cam_id)
        self.max_exhaust = 10  # prevent infinite loop
        self.time_threshold_ms = 10
        self.current_ok = False
        self.current_image = None

    def set(self, *args, **kwargs):
        return self.cap.set(*args, **kwargs)

    def isOpened(self):
        return self.cap.isOpened()

    def read_image(self):
        self.current_ok, self.current_image = self.cap.read()

    def read(self):
        last = None
        for idx in range(self.max_exhaust):
            t = threading.Thread(target=self.read_image)
            t.start()
            t.join(self.time_threshold_ms * 0.001)
            if t.is_alive():
                # still running
                if last is None:
                    # no last frame, wait for this one
                    t.join()
                    return self.current_ok, self.current_image
                else:
                    return last is not None, last
            else:
                # completed
                if self.current_ok:
                    last = self.current_image
                else:
                    last = None
        return False, None


def handle(filename_bytes):
    logging.info("handling ....\n")
    filename = codecs.decode(filename_bytes, 'utf8')
    logging.info("capturing for " + filename + "\n")
    ok, img = capture.read()
    if ok:
        #re.sub(r'[^a-zA-Z]', '_', filename)
        cv2.imwrite(filename, img)
    else:
        logging.info("capture fail\n")


def main():
    global capture

    logging.info("starting with "+" ".join(sys.argv)+"\n")

    error = "Photographer takes four arguments, camera index, width, height and port"
    assert len(sys.argv) == 5, error
    try:
        cam_idx, width, height = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
        port = sys.argv[4]
    except ValueError:
        logging.error(error)
        sys.exit(1)
    assert 0 <= cam_idx < 10, "Camera index is not in the range 0 to 10"

    capture = VCTimeThread(cam_idx)

    logging.info("capture created!")

    capture.set(cv2.CAP_PROP_FRAME_WIDTH, width)
    capture.set(cv2.CAP_PROP_FRAME_HEIGHT, height)

    server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    logging.info("socket created!")
    host = socket.gethostbyname("localhost")
    server_socket.bind((host, int(port)))
    logging.info("bound!")
    server_socket.listen(1)
    # only expect one connection
    logging.info("listening!")

    connection, client_address = server_socket.accept()

    while True:
        data = connection.recv(1024)
        if data:
            handle(data)


if __name__ == "__main__":
    try:
        logging.basicConfig(filename="photo_test.log", level=logging.DEBUG)
        main()
    except:
        logging.error(traceback.format_exc())
        traceback.print_exc()
