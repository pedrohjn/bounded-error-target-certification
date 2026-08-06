function e=sample_disk_error(delta)
r=delta*sqrt(rand); theta=2*pi*rand; e=r*[cos(theta) sin(theta)];
end
