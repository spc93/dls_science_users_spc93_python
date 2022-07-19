#nexus version

def pppscan(scan_no):
    
    n=pdnx(p % scan_no)
    c1n=n.C1/n.ic1monitor;
    c2n=n.C2/n.ic1monitor;
    c3n=n.C3/n.ic1monitor;
    c4n=n.C4/n.ic1monitor;

    plt.figure(); plt.hold(1); 
    plt.plot(n.ppp_offset, n.ic1monitor/max(n.ic1monitor)*2);
    plt.plot(n.ppp_offset, 2*(c1n-min(c1n))/(abs(max(c1n)-min(c1n)))-1,'c'); 
    plt.plot(n.ppp_offset, 2*(c2n-min(c2n))/(abs(max(c2n)-min(c2n)))-1,'r'); 
    plt.plot(n.ppp_offset, 2*(c3n-min(c3n))/(abs(max(c3n)-min(c3n)))-1,'k'); 
    plt.plot(n.ppp_offset, 2*(c4n-min(c4n))/(abs(max(c4n)-min(c4n)))-1,'b'); 
    title(n.scan)
    plt.grid(1); plt.axis('tight')

