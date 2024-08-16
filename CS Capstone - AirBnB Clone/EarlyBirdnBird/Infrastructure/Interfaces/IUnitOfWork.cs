using Infrastructure.Models;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Infrastructure.Interfaces
{
    public interface IUnitOfWork
    {
        public IGenericRepository<Amenity> Amenity { get; }
        public IGenericRepository<Bed> Bed { get; }
        public IGenericRepository<Block> Block { get; }
        public IGenericRepository<Image> Image { get; }
        public IGenericRepository<Price> Price { get; }
        public IGenericRepository<Rental> Rental { get; }
        public IGenericRepository<Reservation> Reservation { get; }
        public IGenericRepository<AppUser> AppUser { get; }
        public IGenericRepository<UserMail> UserMail { get; }
        public IGenericRepository<Review> Review { get; }




        //save changes to the data source

        int Commit();

        Task<int> CommitAsync();
    }
}
