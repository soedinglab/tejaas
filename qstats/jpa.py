class JPA():
    def __init__(self, rank, comm):
        data = None
        if rank == 0:
            data = (1,'a','z',3.14)
            comm.send(data, dest=1, tag=11)
        else:
            print ('on task',rank,'before recv:   data = ',data)
            data = comm.recv(source=0, tag=11)
            print ('on task',rank,'after recv:    data = ',data)
